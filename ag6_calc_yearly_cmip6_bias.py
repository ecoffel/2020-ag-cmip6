import xarray as xr
import xesmf as xe
import numpy as np
import matplotlib.pyplot as plt
import scipy
import statsmodels.api as sm
import cartopy
import cartopy.util
import cartopy.crs as ccrs
import glob
import sys, os
import pickle, gzip
import datetime

dirCmip6 = '/home/edcoffel/drive/MAX-Filer/Research/Climate-02/Data-02-edcoffel-F20/CMIP6'
dirERA5 = '/home/edcoffel/drive/MAX-Filer/Research/Climate-02/Data-02-edcoffel-F20/ERA5'
dirDeepak = '/home/edcoffel/drive/MAX-Filer/Research/Climate-01/Personal-F20/edcoffel-F20/data/projects/ag-land-climate/deepak'
dirAgData = '/home/edcoffel/drive/MAX-Filer/Research/Climate-01/Personal-F20/edcoffel-F20/data/projects/ag-land-climate'
dirProj = '/home/edcoffel/drive/MAX-Filer/Research/Climate-01/Personal-F20/edcoffel-F20/research/2020-ag-cmip6'

with gzip.open('%s/gdd-kdd-lat-era5.dat'%dirAgData, 'rb') as f:
    era5_lat = pickle.load(f)
with gzip.open('%s/gdd-kdd-lon-era5.dat'%dirAgData, 'rb') as f:
    era5_lon = pickle.load(f)
    
region = 'global'
model = sys.argv[1]
recalc = False

if os.path.isfile('cmip6_output/bias/yearly-cmip6-era5-tasmax-max-bias-%s-%s.dat'%(region, model)) and not recalc:
    sys.exit()
    
print('loading regridded tasmax for %s'%model)
cmip6_tasmax_max = xr.open_dataset('cmip6_output/cmip6_tasmax_max_regrid_%s_%s.nc'%(region, model))
cmip6_tasmax_monthly_max = xr.open_dataset('cmip6_output/cmip6_tasmax_monthly_max_regrid_%s_%s.nc'%(region, model))
cmip6_tasmax_mean = xr.open_dataset('cmip6_output/cmip6_tasmax_mean_regrid_%s_%s.nc'%(region, model))


print('loading pre-computed era5...')
era5_tasmax_max_regrid = xr.open_dataset('era5_tasmax_max_regrid_%s.nc'%region)
era5_tasmax_monthly_max_regrid = xr.open_dataset('era5_tasmax_monthly_max_regrid_%s.nc'%region)
era5_tasmax_mean_regrid = xr.open_dataset('era5_tasmax_mean_regrid_%s.nc'%region)

monthly_groups = cmip6_tasmax_monthly_max.groupby('time.month').groups

yearly_tasmax_max_bias = np.full([era5_tasmax_max_regrid.time.values.shape[0], \
                                  era5_tasmax_max_regrid.lat.values.shape[0], \
                                  era5_tasmax_max_regrid.lon.values.shape[0]], np.nan)
yearly_tasmax_monthly_max_bias = np.full([12, len(monthly_groups[1]), \
                                  cmip6_tasmax_monthly_max.lat.values.shape[0], \
                                  cmip6_tasmax_monthly_max.lon.values.shape[0]], np.nan)
yearly_tasmax_mean_bias = np.full([era5_tasmax_mean_regrid.time.values.shape[0], \
                                   era5_tasmax_mean_regrid.lat.values.shape[0], \
                                   era5_tasmax_mean_regrid.lon.values.shape[0]], np.nan)


print('processing %s...'%model)
for y, year in enumerate(range(1981, 2014+1)):
    print('year %d...'%year)
    for xlat in range(yearly_tasmax_max_bias.shape[1]):
        for ylon in range(yearly_tasmax_max_bias.shape[2]):
            yearly_tasmax_max_bias[y, xlat, ylon] = cmip6_tasmax_max.tasmax_max.values[y, xlat, ylon] - \
                                                        era5_tasmax_max_regrid.tasmax_max.values[y, xlat, ylon]
            yearly_tasmax_mean_bias[y, xlat, ylon] = cmip6_tasmax_mean.tasmax_mean.values[y, xlat, ylon] - \
                                                        era5_tasmax_mean_regrid.tasmax_mean.values[y, xlat, ylon]

            for month in range(12):
                yearly_tasmax_monthly_max_bias[month, y, xlat, ylon] = cmip6_tasmax_monthly_max.tasmax_monthly_max.values[monthly_groups[month+1][y], xlat, ylon] - \
                                                        era5_tasmax_monthly_max_regrid.tasmax_monthly_max.values[monthly_groups[month+1][y], xlat, ylon]

with open('cmip6_output/bias/yearly-cmip6-era5-tasmax-max-bias-%s-%s.dat'%(region, model), 'wb') as f:
    pickle.dump(yearly_tasmax_max_bias, f)
with open('cmip6_output/bias/yearly-cmip6-era5-tasmax-monthly-max-bias-%s-%s.dat'%(region, model), 'wb') as f:
    pickle.dump(yearly_tasmax_monthly_max_bias, f)
with open('cmip6_output/bias/yearly-cmip6-era5-tasmax-mean-bias-%s-%s.dat'%(region, model), 'wb') as f:
    pickle.dump(yearly_tasmax_mean_bias, f)