import xarray as xr
import xesmf as xe
import numpy as np
import matplotlib.pyplot as plt
import scipy
import statsmodels.api as sm
import cartopy
import cartopy.crs as ccrs
import glob
import sys
import datetime


gldas_model = 'NOAH'

if gldas_model == 'NOAH':
    dirGLDAS = '/home/edcoffel/drive/MAX-Filer/Research/Climate-01/Data-edcoffel-F20/GLDAS/noah-2-10'
else:
    dirGLDAS = '/home/edcoffel/drive/MAX-Filer/Research/Climate-01/Data-edcoffel-F20/GLDAS/vic-2-10'
    
dirSacks = '/home/edcoffel/drive/MAX-Filer/Research/Climate-01/Personal-F20/edcoffel-F20/data/projects/ag-land-climate'

# transpiration: Tveg_tavg (NOAH) or TVeg_tavg (VIC), soil evap: ESoil_tavg, canopy evap: ECanop_tavg, et: Evap_tavg 
orig_var = 'Tveg_tavg'
crop = 'Maize'

year = int(sys.argv[1])

latRange = [-90, 90]
lonRange = [0, 360]


# LOADING SACKS CROP CALENDARS
sacksMaizeNc = xr.open_dataset('%s/sacks/%s.crop.calendar.fill.nc'%(dirSacks, crop))
sacksStart = sacksMaizeNc['plant'].values
sacksStart = np.roll(sacksStart, -int(sacksStart.shape[1]/2), axis=1)
sacksStart[sacksStart < 0] = np.nan
sacksEnd = sacksMaizeNc['harvest'].values
sacksEnd = np.roll(sacksEnd, -int(sacksEnd.shape[1]/2), axis=1)
sacksEnd[sacksEnd < 0] = np.nan

# THESE ARE THE LAT/LON GRIDS THAT SACKS IS ON
sacksLat = np.linspace(90, -90, 360)
sacksLon = np.linspace(0, 360, 720)

regridMesh_cur_model = xr.Dataset()

n = 0

# LOAD 1 YEAR OF ERA5 DATA
print('opening gldas %d...'%year)
gldas = xr.open_mfdataset('%s/GLDAS_%s10_M.A%d*.020.nc4'%(dirGLDAS, gldas_model, year))
gldas.load()

print('opening gldas %d...'%(year-1))
gldas_last_year = xr.open_mfdataset('%s/GLDAS_%s10_M.A%d*.020.nc4'%(dirGLDAS, gldas_model, year-1))
gldas_last_year.load()

if orig_var != 'Evap_tavg':
    # W/m2 -> kg/(m2*s) = mm/s
    gldas[orig_var] /= (2260*1000)
    gldas_last_year[orig_var] /= (2260*1000)
    # mm/s -> mm/day
    gldas[orig_var] *= 60*60*24
    gldas_last_year[orig_var] *= 60*60*24
else:
    # mm/s -> mm/day
    gldas[orig_var] *= 60*60*24
    gldas_last_year[orig_var] *= 60*60*24

# THIS USES XESMF TO REGRID THE SACKS DATA TO ERA5 RES
# regrid sacks data to current model res
regridMesh_cur_model = xr.Dataset({'lat': (['lat'], sacksLat),
                                   'lon': (['lon'], sacksLon)})

regridder_gldas = xe.Regridder(xr.DataArray(data=gldas[orig_var], dims=['time', 'lat', 'lon'], coords={'lat':gldas.lat, 'lon':gldas.lon}), regridMesh_cur_model, 'bilinear', reuse_weights=True)
regridder_gldas_last_year = xe.Regridder(xr.DataArray(data=gldas_last_year[orig_var], dims=['time', 'lat', 'lon'], coords={'lat':gldas_last_year.lat, 'lon':gldas_last_year.lon}), regridMesh_cur_model, 'bilinear', reuse_weights=True)

gldas = regridder_gldas(gldas[orig_var])
gldas_last_year = regridder_gldas_last_year(gldas_last_year[orig_var])


# count up all non-nan grid cells so we can estimate percent complete
ngrid = 0
for xlat in range(gldas.lat.size):
    for ylon in range(gldas.lon.size):
        
        if not np.isnan(sacksStart[xlat, ylon]):
            curStart = datetime.datetime.strptime('2020%d'%(round(sacksStart[xlat, ylon])+1), '%Y%j').date().month
            sacksStart[xlat, ylon] = curStart-1
            
        if not np.isnan(sacksEnd[xlat, ylon]):
            curEnd = datetime.datetime.strptime('2020%d'%(round(sacksEnd[xlat, ylon])+1), '%Y%j').date().month
            sacksEnd[xlat, ylon] = curEnd-1
        
        if ~np.isnan(sacksStart[xlat, ylon]) and ~np.isnan(sacksEnd[xlat, ylon]):
            ngrid += 1


yearly_grow_mean = np.full([gldas.lat.size, gldas.lon.size], np.nan)
            
            
# THIS LOOPS OVER EVERY GRID CELL OF ERA5 AND EXTRACTS DAILY ERA5 DATA THAT FALLS WITHIN THE SACKS GROWING SEASON
# latitude loop
for xlat in range(gldas.lat.size):
    # longitude loop
    for ylon in range(gldas.lon.size):

        # if sacks calendar is defined at this grid cell
        if ~np.isnan(sacksStart[xlat, ylon]) and ~np.isnan(sacksEnd[xlat, ylon]):
    
            # just print out our progress
            if n % 1000 == 0:
                print('%.2f%%'%(n/ngrid*100))

            # there are 2 possibilities - that the planting date is before the harvest date in the current year (northern hemisphere), 
            # or that the planting date is late in the year and the harvest date is in the beginning of the next year (southern hemisphere)
            # we need to handle these two cases separately
            
            # SOUTHERN HEMISPHERE - NEED 2 YEARS OF ERA5 DATA (this year and last year)
            if sacksStart[xlat, ylon] > sacksEnd[xlat, ylon]:

                # start loop on 2nd year to allow for growing season that crosses jan 1
                curVar1 = gldas_last_year[int(sacksStart[xlat, ylon]):, xlat, ylon].values
                curVar2 = gldas[:int(sacksEnd[xlat, ylon]), xlat, ylon].values

                curVar = np.concatenate([curVar1, curVar2])

                if len(curVar) > 0:
                    yearly_grow_mean[xlat, ylon] = np.nanmean(curVar)
                n += 1
                

            # NORTHERN HEMISPHERE - SIMPLER, JUST NEED 1 YEAR OF ERA5
            else:
                curVar = gldas[int(sacksStart[xlat, ylon]):int(sacksEnd[xlat, ylon]), xlat, ylon].values
                if len(curVar) > 0:
                    yearly_grow_mean[xlat, ylon] = np.nanmean(curVar)
                    
                n += 1

print('renaming dims...')

# SAVE THE EXTRACTED TEMP DATA
da_grow_mean = xr.DataArray(data   = yearly_grow_mean, 
                      dims   = ['lat', 'lon'],
                      coords = {'time': year, 'lat':gldas.lat, 'lon':gldas.lon},
                      attrs  = {'units'     : 'mm'
                        })
ds_grow_mean = xr.Dataset()
ds_grow_mean['evap_grow_mean'] = da_grow_mean

print('saving netcdf...')
ds_grow_mean.to_netcdf('gldas_output/gldas-%s_%s_%s_%s_grow_mean_global_%d_fixed_sh.nc'%(gldas_model, orig_var, crop, orig_var, year))
