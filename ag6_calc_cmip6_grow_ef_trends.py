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
dirProj = '/home/edcoffel/drive/MAX-Filer/Research/Climate-01/Personal-F20/edcoffel-F20/research/2020-ag-cmip6'

model = sys.argv[1]
region = 'global'
crop = 'Maize'

if not os.path.isfile('cmip6_output/growing_season/cmip6_%s_grow_ef_mon_%s_%s_regrid.nc'%(crop, region, model)):
    print('no base file, quitting')
    sys.exit()

cmip6_ef = xr.open_dataset('cmip6_output/growing_season/cmip6_%s_grow_ef_mon_%s_%s_regrid.nc'%(crop, region, model))

print('calc trends for %s...'%model)
tmp_ef_trend = np.full([cmip6_ef.lat.size, cmip6_ef.lon.size], np.nan)

for xlat in range(tmp_ef_trend.shape[0]):
    for ylon in range(tmp_ef_trend.shape[1]):
        
        
        cur_ef = cmip6_ef['grow_ef'].values[:, xlat, ylon]
        nn = np.where(~np.isnan(cur_ef))[0]
        if len(nn) > 10:
            X = sm.add_constant(range(len(nn)))
            mdl = sm.OLS(cur_ef[nn], X).fit()
            tmp_ef_trend[xlat, ylon] = mdl.params[1]*10

tmpDs_ef = xr.DataArray(data   = tmp_ef_trend, 
                  dims   = ['lat', 'lon'],
                  coords = {'lat':cmip6_ef.lat, 'lon':cmip6_ef.lon},
                  attrs  = {'units'     : 'Fraction'
                    })


cmip6_ef_trend = xr.Dataset()

cmip6_ef_trend['grow_ef'] = tmpDs_ef.assign_coords({'model':model})

cmip6_ef_trend.to_netcdf('cmip6_output/growing_season/cmip6_%s_grow_ef_mon_trend_%s_%s_regrid.nc'%(crop, region, model))
