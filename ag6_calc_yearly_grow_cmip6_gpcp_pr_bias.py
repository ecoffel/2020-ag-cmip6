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

print('loading regridded pr for %s'%model)
cmip6_pr_grow = xr.open_dataset('cmip6_output/growing_season/cmip6_%s_grow_pr_mean_%s_%s_regrid_fixed_sh.nc'%(crop, region, model))

print('loading pre-computed era5...')
gpcp_pr_grow_regrid = xr.open_dataset('gpcp_output/gpcp_Maize_grow_mean_global_regrid_fixed_sh.nc')


yearly_pr_grow_bias = np.full([len(range(1981, 2014+1)), \
                                   gpcp_pr_grow_regrid.lat.values.shape[0], \
                                   gpcp_pr_grow_regrid.lon.values.shape[0]], np.nan)


print('processing %s...'%model)
for y, year in enumerate(range(1981, 2014+1)):
    print('year %d...'%year)
    for xlat in range(yearly_pr_grow_bias.shape[1]):
        for ylon in range(yearly_pr_grow_bias.shape[2]):
            yearly_pr_grow_bias[y, xlat, ylon] = cmip6_pr_grow.pr_grow_mean.values[y, xlat, ylon] - \
                                                        gpcp_pr_grow_regrid.precip_grow_mean.values[y, xlat, ylon]

with open('cmip6_output/bias/yearly-cmip6-gpcp-pr-grow-bias-%s-%s-fixed-sh.dat'%(region, model), 'wb') as f:
    pickle.dump(yearly_pr_grow_bias, f)