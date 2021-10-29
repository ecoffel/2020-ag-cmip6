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
dirAgData = '/home/edcoffel/drive/MAX-Filer/Research/Climate-01/Personal-F20/edcoffel-F20/data/projects/ag-land-climate'
dirProj = '/home/edcoffel/drive/MAX-Filer/Research/Climate-01/Personal-F20/edcoffel-F20/research/2020-ag-cmip6'

with gzip.open('%s/gdd-kdd-lat-era5.dat'%dirAgData, 'rb') as f:
    era5_lat = pickle.load(f)
with gzip.open('%s/gdd-kdd-lon-era5.dat'%dirAgData, 'rb') as f:
    era5_lon = pickle.load(f)
    
crop = 'Maize'
region = 'global'
model = sys.argv[1]

cmip6_var = 'lai'
era5_var = 'lai'

lai_obs_regrid = xr.open_dataset('lai_data/lai_%s_grow_mean_regrid_global.nc'%(crop))

print('loading regridded lai for %s'%model)
cmip6_lai_grow = xr.open_dataset('cmip6_output/growing_season/cmip6_%s_grow_%s_mon_global_%s_regrid.nc'%(crop, cmip6_var, model))


lai_grow_bias = np.full([lai_obs_regrid.lat.values.shape[0], \
                                   lai_obs_regrid.lon.values.shape[0]], np.nan)


print('processing %s...'%model)
for xlat in range(lai_grow_bias.shape[0]):
    for ylon in range(lai_grow_bias.shape[1]):
        lai_grow_bias[xlat, ylon] = np.nanmean(cmip6_lai_grow['%s_grow_mean'%cmip6_var].values[:, xlat, ylon]) - \
                                                    lai_obs_regrid.lai_grow_mean.values[xlat, ylon]

with open('cmip6_output/bias/climo-cmip6-obs-lai-grow-bias-%s-%s.dat'%(region, model), 'wb') as f:
    pickle.dump(lai_grow_bias, f)