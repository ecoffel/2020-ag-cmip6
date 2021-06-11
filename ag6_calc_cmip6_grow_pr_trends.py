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

if not os.path.isfile('cmip6_output/growing_season/cmip6_%s_grow_pr_mean_%s_%s_regrid.nc'%(crop, region, model)):
    print('no base file, quitting')
    sys.exit()

cmip6_pr = xr.open_dataset('cmip6_output/growing_season/cmip6_%s_grow_pr_mean_%s_%s_regrid.nc'%(crop, region, model))

print('calc trends for %s...'%model)
tmp_pr_trend = np.full([cmip6_pr.lat.size, cmip6_pr.lon.size], np.nan)

for xlat in range(tmp_pr_trend.shape[0]):
    for ylon in range(tmp_pr_trend.shape[1]):
        
        
        cur_pr = cmip6_pr['pr_grow_mean'].values[:, xlat, ylon]
        nn = np.where(~np.isnan(cur_pr))[0]
        if len(nn) > 10:
            X = sm.add_constant(range(len(nn)))
            mdl = sm.OLS(cur_pr[nn], X).fit()
            tmp_pr_trend[xlat, ylon] = mdl.params[1]*10

tmpDs_pr = xr.DataArray(data   = tmp_pr_trend, 
                  dims   = ['lat', 'lon'],
                  coords = {'lat':cmip6_pr.lat, 'lon':cmip6_pr.lon},
                  attrs  = {'units'     : 'mm/yr/decade'
                    })


cmip6_pr_trend = xr.Dataset()

cmip6_pr_trend['pr_grow_mean'] = tmpDs_pr.assign_coords({'model':model})

cmip6_pr_trend.to_netcdf('cmip6_output/growing_season/cmip6_%s_grow_pr_mean_trend_%s_%s_regrid.nc'%(crop, region, model))
