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
    
crop = 'Maize'
region = 'global'
rcp = 'historical'
model = sys.argv[1]
member = sys.argv[2]

print('loading regridded tasmax for %s'%model)
cmip6_tasmax_grow_max = xr.open_dataset('cmip6_output/growing_season/cmip6_%s_tasmax_grow_max_%s_%s_%s_regrid.nc'%(crop, region, model, member))


print('loading pre-computed era5...')
era5_tasmax_grow_max_regrid = xr.open_dataset('era5/growing_season/era5_%s_tasmax_grow_max_regrid_%s.nc'%(crop,region))

yearly_tasmax_grow_max_bias = np.full([len(range(1981, 2014+1)), \
                                  era5_tasmax_grow_max_regrid.lat.values.shape[0], \
                                  era5_tasmax_grow_max_regrid.lon.values.shape[0]], np.nan)



print('processing %s...'%model)
for y, year in enumerate(range(1981, 2014+1)):
    print('year %d...'%year)
    for xlat in range(yearly_tasmax_grow_max_bias.shape[1]):
        for ylon in range(yearly_tasmax_grow_max_bias.shape[2]):
            yearly_tasmax_grow_max_bias[y, xlat, ylon] = cmip6_tasmax_grow_max.tasmax_grow_max.values[y, xlat, ylon] - \
                                                        era5_tasmax_grow_max_regrid.tasmax_grow_max.values[y, xlat, ylon]

with open('cmip6_output/bias/yearly-cmip6-era5-tasmax-grow-max-bias-%s-%s-%s.dat'%(region, model, member), 'wb') as f:
    pickle.dump(yearly_tasmax_grow_max_bias, f)
