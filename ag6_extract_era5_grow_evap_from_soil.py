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

dirERA5Land = '/home/edcoffel/drive/MAX-Filer/Research/Climate-02/Data-02-edcoffel-F20/ERA5-Land'
dirSacks = '/home/edcoffel/drive/MAX-Filer/Research/Climate-01/Personal-F20/edcoffel-F20/data/projects/ag-land-climate'

file_var = 'evap_from_bare_soil'
orig_var = 'evaow' # this is the variable for evap from water because the era5 vars are switched: https://confluence.ecmwf.int/display/CKB/ERA5-Land%3A+data+documentation
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
print('opening era5 %d...'%year)
evap_era5 = xr.open_dataset('%s/monthly/%s_%d.nc'%(dirERA5Land, file_var, year))
evap_era5.load()

print('opening era5 %d...'%(year-1))
evap_era5_last_year = xr.open_dataset('%s/monthly/%s_%d.nc'%(dirERA5Land, file_var, year-1))
evap_era5_last_year.load()

evap_era5 = evap_era5.rename_dims(latitude='lat', longitude='lon')
evap_era5 = evap_era5.rename({'latitude':'lat', 'longitude':'lon'})

evap_era5_last_year = evap_era5_last_year.rename_dims(latitude='lat', longitude='lon')
evap_era5_last_year = evap_era5_last_year.rename({'latitude':'lat', 'longitude':'lon'})

# THIS USES XESMF TO REGRID THE SACKS DATA TO ERA5 RES
# regrid sacks data to current model res
regridMesh_cur_model = xr.Dataset({'lat': (['lat'], sacksLat),
                                   'lon': (['lon'], sacksLon)})

regridder_evap = xe.Regridder(xr.DataArray(data=evap_era5[orig_var], dims=['time', 'lat', 'lon'], coords={'lat':evap_era5.lat, 'lon':evap_era5.lon}), regridMesh_cur_model, 'bilinear', reuse_weights=True)
regridder_evap_last_year = xe.Regridder(xr.DataArray(data=evap_era5_last_year[orig_var], dims=['time', 'lat', 'lon'], coords={'lat':evap_era5_last_year.lat, 'lon':evap_era5_last_year.lon}), regridMesh_cur_model, 'bilinear', reuse_weights=True)

evap_era5 = regridder_evap(evap_era5)
evap_era5_last_year = regridder_evap_last_year(evap_era5_last_year)

# count up all non-nan grid cells so we can estimate percent complete
ngrid = 0
for xlat in range(evap_era5.lat.size):
    for ylon in range(evap_era5.lon.size):
        
        if not np.isnan(sacksStart[xlat, ylon]):
            curStart = datetime.datetime.strptime('2020%d'%(round(sacksStart[xlat, ylon])+1), '%Y%j').date().month
            sacksStart[xlat, ylon] = curStart
            
        if not np.isnan(sacksEnd[xlat, ylon]):
            curEnd = datetime.datetime.strptime('2020%d'%(round(sacksEnd[xlat, ylon])+1), '%Y%j').date().month
            sacksEnd[xlat, ylon] = curEnd
        
        if ~np.isnan(sacksStart[xlat, ylon]) and ~np.isnan(sacksEnd[xlat, ylon]):
            ngrid += 1


yearly_grow_evap_mean = np.full([evap_era5.lat.size, evap_era5.lon.size], np.nan)
            
            
# THIS LOOPS OVER EVERY GRID CELL OF ERA5 AND EXTRACTS DAILY ERA5 DATA THAT FALLS WITHIN THE SACKS GROWING SEASON
# latitude loop
for xlat in range(evap_era5.lat.size):
    # longitude loop
    for ylon in range(evap_era5.lon.size):

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
                curevap1 = evap_era5_last_year[orig_var][int(sacksEnd[xlat, ylon]):, xlat, ylon]
                curevap2 = evap_era5[orig_var][:int(sacksStart[xlat, ylon]), xlat, ylon]

                curevap = np.concatenate([curevap1, curevap2])

                if len(curevap) > 0:
                    yearly_grow_evap_mean[xlat, ylon] = np.nanmean(curevap)
                n += 1
                

            # NORTHERN HEMISPHERE - SIMPLER, JUST NEED 1 YEAR OF ERA5
            else:
                curevap = evap_era5[orig_var][int(sacksStart[xlat, ylon]):int(sacksEnd[xlat, ylon]), xlat, ylon]
                if len(curevap) > 0:
                    yearly_grow_evap_mean[xlat, ylon] = np.nanmean(curevap)
                    
                n += 1

print('renaming dims...')

# SAVE THE EXTRACTED TEMP DATA
da_grow_evap_mean = xr.DataArray(data   = yearly_grow_evap_mean, 
                      dims   = ['lat', 'lon'],
                      coords = {'time': year, 'lat':evap_era5.lat, 'lon':evap_era5.lon},
                      attrs  = {'units'     : 'mm'
                        })
ds_grow_evap_mean = xr.Dataset()
ds_grow_evap_mean['evap_grow_mean'] = da_grow_evap_mean

print('saving netcdf...')
ds_grow_evap_mean.to_netcdf('era5/growing_season/era5_%s_evap_from_soil_grow_mean_global_%d.nc'%(crop, year))
