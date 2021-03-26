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

file_var = 'tasmax'
orig_var = 'mx2t'

latRange = [-90, 90]
lonRange = [0, 360]

print('opening era5...')
temp_era5 = xr.open_mfdataset('%s/daily/%s_*.nc'%(dirERA5, file_var))

print('resampling era5...')
temp_max_era5 = temp_era5.resample(time='1Y').max(dim='time')
temp_max_era5[orig_var] -= 273.15

temp_median_era5 = temp_era5.resample(time='1Y').median(dim='time')
temp_median_era5[orig_var] -= 273.15

print('renaming dims...')
temp_max_era5 = temp_max_era5.rename_dims(latitude='lat', longitude='lon')
temp_median_era5 = temp_median_era5.rename_dims(latitude='lat', longitude='lon')

temp_max_era5 = temp_max_era5.rename({'latitude':'lat', 'longitude':'lon', orig_var:'%s_max'%file_var})
temp_median_era5 = temp_median_era5.rename({'latitude':'lat', 'longitude':'lon', orig_var:'%s_median'%file_var})

print('writing era5 txx...')
temp_max_era5.to_netcdf('era5_%s_max_global.nc'%file_var)
print('writing era5 t50p...')
temp_median_era5.to_netcdf('era5_%s_median_global.nc'%file_var)