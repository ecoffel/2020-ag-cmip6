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

dirERA5 = '/home/edcoffel/drive/MAX-Filer/Research/Climate-02/Data-02-edcoffel-F20/ERA5'
dirSacks = '/home/edcoffel/drive/MAX-Filer/Research/Climate-01/Personal-F20/edcoffel-F20/data/projects/ag-land-climate'

file_var = 'tp'
orig_var = 'tp'
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
pr_era5 = xr.open_dataset('%s/daily/%s_%d.nc'%(dirERA5, file_var, year))
pr_era5.load()

print('opening era5 %d...'%(year-1))
pr_era5_last_year = xr.open_dataset('%s/daily/%s_%d.nc'%(dirERA5, file_var, year-1))
pr_era5_last_year.load()

pr_era5 = pr_era5.rename_dims(latitude='lat', longitude='lon')
pr_era5 = pr_era5.rename({'latitude':'lat', 'longitude':'lon'})

pr_era5_last_year = pr_era5_last_year.rename_dims(latitude='lat', longitude='lon')
pr_era5_last_year = pr_era5_last_year.rename({'latitude':'lat', 'longitude':'lon'})

# THIS USES XESMF TO REGRID THE SACKS DATA TO ERA5 RES
# regrid sacks data to current model res
regridMesh_cur_model = xr.Dataset({'lat': (['lat'], pr_era5.lat),
                                   'lon': (['lon'], pr_era5.lon)})

regridder_start = xe.Regridder(xr.DataArray(data=sacksStart, dims=['lat', 'lon'], coords={'lat':sacksLat, 'lon':sacksLon}), regridMesh_cur_model, 'bilinear', reuse_weights=True)
regridder_end = xe.Regridder(xr.DataArray(data=sacksEnd, dims=['lat', 'lon'], coords={'lat':sacksLat, 'lon':sacksLon}), regridMesh_cur_model, 'bilinear', reuse_weights=True)

sacksStart_regrid = regridder_start(sacksStart)
sacksEnd_regrid = regridder_end(sacksEnd)

# count up all non-nan grid cells so we can estimate percent complete
ngrid = 0
for xlat in range(pr_era5.lat.size):
    for ylon in range(pr_era5.lon.size):
        if ~np.isnan(sacksStart_regrid[xlat, ylon]) and ~np.isnan(sacksEnd_regrid[xlat, ylon]):
            ngrid += 1


yearly_grow_pr_mean = np.full([pr_era5.lat.size, pr_era5.lon.size], np.nan)
            
            
# THIS LOOPS OVER EVERY GRID CELL OF ERA5 AND EXTRACTS DAILY ERA5 DATA THAT FALLS WITHIN THE SACKS GROWING SEASON
# latitude loop
for xlat in range(pr_era5.lat.size):
    # longitude loop
    for ylon in range(pr_era5.lon.size):

        # if sacks calendar is defined at this grid cell
        if ~np.isnan(sacksStart_regrid[xlat, ylon]) and ~np.isnan(sacksEnd_regrid[xlat, ylon]):
    
            # just print out our progress
            if n % 1000 == 0:
                print('%.2f%%'%(n/ngrid*100))

            # there are 2 possibilities - that the planting date is before the harvest date in the current year (northern hemisphere), 
            # or that the planting date is late in the year and the harvest date is in the beginning of the next year (southern hemisphere)
            # we need to handle these two cases separately
            
            # SOUTHERN HEMISPHERE - NEED 2 YEARS OF ERA5 DATA (this year and last year)
            if sacksStart_regrid[xlat, ylon] > sacksEnd_regrid[xlat, ylon]:

                # start loop on 2nd year to allow for growing season that crosses jan 1
                curPr1 = pr_era5_last_year[orig_var][int(sacksEnd_regrid[xlat, ylon]):, xlat, ylon]
                curPr2 = pr_era5[orig_var][:int(sacksStart_regrid[xlat, ylon]), xlat, ylon]

                curPr = np.concatenate([curPr1, curPr2])

                if len(curPr) > 0:
                    yearly_grow_pr_mean[xlat, ylon] = np.nanmean(curPr)
                n += 1
                

            # NORTHERN HEMISPHERE - SIMPLER, JUST NEED 1 YEAR OF ERA5
            else:
                curPr = pr_era5[orig_var][int(sacksStart_regrid[xlat, ylon]):int(sacksEnd_regrid[xlat, ylon]), xlat, ylon]
                if len(curPr) > 0:
                    yearly_grow_pr_mean[xlat, ylon] = np.nanmean(curPr)
                n += 1

print('renaming dims...')

# SAVE THE EXTRACTED TEMP DATA
da_grow_pr_mean = xr.DataArray(data   = yearly_grow_pr_mean, 
                      dims   = ['lat', 'lon'],
                      coords = {'time': year, 'lat':pr_era5.lat, 'lon':pr_era5.lon},
                      attrs  = {'units'     : 'mm'
                        })
ds_grow_pr_mean = xr.Dataset()
ds_grow_pr_mean['pr_grow_mean'] = da_grow_pr_mean

print('saving netcdf...')
ds_grow_pr_mean.to_netcdf('era5/growing_season/era5_%s_pr_grow_mean_global_%d.nc'%(crop, year))
