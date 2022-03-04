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
dirERA5 = '/home/edcoffel/drive/MAX-Filer/Research/Climate-02/Data-02-edcoffel-F20/ERA5'
dirSacks = '/home/edcoffel/drive/MAX-Filer/Research/Climate-01/Personal-F20/edcoffel-F20/data/projects/ag-land-climate'

file_var1 = 'lai_low'
file_var2 = 'lai_high'
orig_var1 = 'lai_lv'
orig_var2 = 'lai_hv'
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

# regrid sacks data to current model res
regridMesh_cur_model = xr.Dataset({'lat': (['lat'], sacksLat),
                                   'lon': (['lon'], sacksLon)})

n = 0

# load low vegetaion cover
era5_low_veg = xr.open_dataset('%s/monthly/low_vegetation_cover.nc'%dirERA5)

era5_low_veg = era5_low_veg.rename_dims(latitude='lat', longitude='lon')
era5_low_veg = era5_low_veg.rename({'latitude':'lat', 'longitude':'lon'})

era5_low_veg = era5_low_veg.isel(time=0,expver=0)
era5_low_veg.load()

era5_low_veg = era5_low_veg.drop('expver')
era5_low_veg = era5_low_veg.drop('time')

regridder_low_veg = xe.Regridder(xr.DataArray(data=era5_low_veg.cvl, dims=['lat', 'lon'], coords={'lat':era5_low_veg.lat, 'lon':era5_low_veg.lon}), regridMesh_cur_model, 'bilinear', reuse_weights=True)
era5_low_veg = regridder_low_veg(era5_low_veg)

# and load high veg cover
era5_high_veg = xr.open_dataset('%s/monthly/high_vegetation_cover.nc'%dirERA5)

era5_high_veg = era5_high_veg.rename_dims(latitude='lat', longitude='lon')
era5_high_veg = era5_high_veg.rename({'latitude':'lat', 'longitude':'lon'})

era5_high_veg = era5_high_veg.isel(time=0,expver=0)
era5_high_veg.load()

era5_high_veg = era5_high_veg.drop('expver')
era5_high_veg = era5_high_veg.drop('time')

regridder_high_veg = xe.Regridder(xr.DataArray(data=era5_high_veg.cvh, dims=['lat', 'lon'], coords={'lat':era5_high_veg.lat, 'lon':era5_high_veg.lon}), regridMesh_cur_model, 'bilinear', reuse_weights=True)
era5_high_veg = regridder_high_veg(era5_high_veg)


# LOAD 1 YEAR OF ERA5 DATA
print('opening era5 %d...'%year)
lai_low_era5 = xr.open_dataset('%s/monthly/%s_%d.nc'%(dirERA5Land, file_var1, year))
lai_low_era5.load()

lai_high_era5 = xr.open_dataset('%s/monthly/%s_%d.nc'%(dirERA5Land, file_var2, year))
lai_high_era5.load()

print('opening era5 %d...'%(year-1))
lai_low_era5_last_year = xr.open_dataset('%s/monthly/%s_%d.nc'%(dirERA5Land, file_var1, year-1))
lai_low_era5_last_year.load()

lai_high_era5_last_year = xr.open_dataset('%s/monthly/%s_%d.nc'%(dirERA5Land, file_var2, year-1))
lai_high_era5_last_year.load()

lai_low_era5 = lai_low_era5.rename_dims(latitude='lat', longitude='lon')
lai_low_era5 = lai_low_era5.rename({'latitude':'lat', 'longitude':'lon'})

lai_high_era5 = lai_high_era5.rename_dims(latitude='lat', longitude='lon')
lai_high_era5 = lai_high_era5.rename({'latitude':'lat', 'longitude':'lon'})

lai_low_era5_last_year = lai_low_era5_last_year.rename_dims(latitude='lat', longitude='lon')
lai_low_era5_last_year = lai_low_era5_last_year.rename({'latitude':'lat', 'longitude':'lon'})

lai_high_era5_last_year = lai_high_era5_last_year.rename_dims(latitude='lat', longitude='lon')
lai_high_era5_last_year = lai_high_era5_last_year.rename({'latitude':'lat', 'longitude':'lon'})

# THIS USES XESMF TO REGRID THE SACKS DATA TO ERA5 RES


regridder_lai_low = xe.Regridder(xr.DataArray(data=lai_low_era5[orig_var1], dims=['time', 'lat', 'lon'], coords={'lat':lai_low_era5.lat, 'lon':lai_low_era5.lon}), regridMesh_cur_model, 'bilinear', reuse_weights=True)
regridder_lai_high = xe.Regridder(xr.DataArray(data=lai_high_era5[orig_var2], dims=['time', 'lat', 'lon'], coords={'lat':lai_high_era5.lat, 'lon':lai_high_era5.lon}), regridMesh_cur_model, 'bilinear', reuse_weights=True)
regridder_lai_low_last_year = xe.Regridder(xr.DataArray(data=lai_low_era5_last_year[orig_var1], dims=['time', 'lat', 'lon'], coords={'lat':lai_low_era5_last_year.lat, 'lon':lai_low_era5_last_year.lon}), \
                                           regridMesh_cur_model, 'bilinear', reuse_weights=True)
regridder_lai_high_last_year = xe.Regridder(xr.DataArray(data=lai_high_era5_last_year[orig_var2], dims=['time', 'lat', 'lon'], coords={'lat':lai_high_era5_last_year.lat, 'lon':lai_high_era5_last_year.lon}), \
                                            regridMesh_cur_model, 'bilinear', reuse_weights=True)

lai_low_era5 = regridder_lai_low(lai_low_era5)
lai_high_era5 = regridder_lai_high(lai_high_era5)
lai_low_era5_last_year = regridder_lai_low_last_year(lai_low_era5_last_year)
lai_high_era5_last_year = regridder_lai_high_last_year(lai_high_era5_last_year)

# count up all non-nan grid cells so we can estimate percent complete
ngrid = 0
for xlat in range(lai_low_era5.lat.size):
    for ylon in range(lai_low_era5.lon.size):
        
        if not np.isnan(sacksStart[xlat, ylon]):
            curStart = datetime.datetime.strptime('2020%d'%(round(sacksStart[xlat, ylon])+1), '%Y%j').date().month
            sacksStart[xlat, ylon] = curStart-1
            
        if not np.isnan(sacksEnd[xlat, ylon]):
            curEnd = datetime.datetime.strptime('2020%d'%(round(sacksEnd[xlat, ylon])+1), '%Y%j').date().month
            sacksEnd[xlat, ylon] = curEnd-1
        
        if ~np.isnan(sacksStart[xlat, ylon]) and ~np.isnan(sacksEnd[xlat, ylon]):
            ngrid += 1


yearly_grow_lai_low_mean = np.full([lai_low_era5.lat.size, lai_low_era5.lon.size], np.nan)
yearly_grow_lai_high_mean = np.full([lai_high_era5.lat.size, lai_high_era5.lon.size], np.nan)
            
            
# THIS LOOPS OVER EVERY GRID CELL OF ERA5 AND EXTRACTS DAILY ERA5 DATA THAT FALLS WITHIN THE SACKS GROWING SEASON
# latitude loop
for xlat in range(lai_low_era5.lat.size):
    # longitude loop
    for ylon in range(lai_low_era5.lon.size):

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
                cur_lai_low1 = lai_low_era5_last_year[orig_var1][int(sacksStart[xlat, ylon]):, xlat, ylon]
                cur_lai_low2 = lai_low_era5[orig_var1][:int(sacksEnd[xlat, ylon]), xlat, ylon]
                
                cur_lai_high1 = lai_high_era5_last_year[orig_var2][int(sacksStart[xlat, ylon]):, xlat, ylon]
                cur_lai_high2 = lai_high_era5[orig_var2][:int(sacksEnd[xlat, ylon]), xlat, ylon]

                cur_lai_low = np.concatenate([cur_lai_low1, cur_lai_low2]) * era5_low_veg.cvl.values[xlat, ylon]
                cur_lai_high = np.concatenate([cur_lai_high1, cur_lai_high2]) * era5_high_veg.cvh.values[xlat, ylon]

                if len(cur_lai_low) > 0 and len(cur_lai_high) > 0:
                    yearly_grow_lai_low_mean[xlat, ylon] = np.nanmean(cur_lai_low)
                    yearly_grow_lai_high_mean[xlat, ylon] = np.nanmean(cur_lai_high)
                n += 1
                

            # NORTHERN HEMISPHERE - SIMPLER, JUST NEED 1 YEAR OF ERA5
            else:
                cur_lai_low = lai_low_era5[orig_var1][int(sacksStart[xlat, ylon]):int(sacksEnd[xlat, ylon]), xlat, ylon]
                cur_lai_high = lai_high_era5[orig_var2][int(sacksStart[xlat, ylon]):int(sacksEnd[xlat, ylon]), xlat, ylon]
                
                if len(cur_lai_low) > 0 and len(cur_lai_high) > 0:
                    yearly_grow_lai_low_mean[xlat, ylon] = np.nanmean(cur_lai_low) * era5_low_veg.cvl.values[xlat, ylon]
                    yearly_grow_lai_high_mean[xlat, ylon] = np.nanmean(cur_lai_high) * era5_high_veg.cvh.values[xlat, ylon]
                    
                n += 1

print('renaming dims...')

# SAVE THE EXTRACTED TEMP DATA
da_grow_lai_mean = xr.DataArray(data   = yearly_grow_lai_low_mean + yearly_grow_lai_high_mean, 
                      dims   = ['lat', 'lon'],
                      coords = {'time': year, 'lat':lai_low_era5.lat, 'lon':lai_low_era5.lon},
                      attrs  = {'units'     : 'LAI'
                        })
ds_grow_lai_mean = xr.Dataset()
ds_grow_lai_mean['lai_grow_mean'] = da_grow_lai_mean

print('saving netcdf...')
ds_grow_lai_mean.to_netcdf('era5/growing_season/era5_%s_lai_grow_mean_global_%d_fixed_sh.nc'%(crop, year))



da_grow_lai_low = xr.DataArray(data   = yearly_grow_lai_low_mean, 
                      dims   = ['lat', 'lon'],
                      coords = {'time': year, 'lat':lai_low_era5.lat, 'lon':lai_low_era5.lon},
                      attrs  = {'units'     : 'LAI'
                        })
ds_grow_lai_low = xr.Dataset()
ds_grow_lai_low['lai_grow_mean'] = da_grow_lai_low

print('saving netcdf...')
ds_grow_lai_low.to_netcdf('era5/growing_season/era5_%s_lai_low_grow_mean_global_%d_fixed_sh.nc'%(crop, year))




da_grow_lai_high = xr.DataArray(data   = yearly_grow_lai_high_mean, 
                      dims   = ['lat', 'lon'],
                      coords = {'time': year, 'lat':lai_high_era5.lat, 'lon':lai_high_era5.lon},
                      attrs  = {'units'     : 'LAI'
                        })
ds_grow_lai_high = xr.Dataset()
ds_grow_lai_high['lai_grow_mean'] = da_grow_lai_high

print('saving netcdf...')
ds_grow_lai_high.to_netcdf('era5/growing_season/era5_%s_lai_high_grow_mean_global_%d_fixed_sh.nc'%(crop, year))
