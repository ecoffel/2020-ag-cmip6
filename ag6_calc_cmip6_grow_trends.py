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
var = 'ef'

if not os.path.isfile('cmip6_output/growing_season/cmip6_%s_tasmax_grow_max_%s_%s_regrid.nc'%(crop, region, model)):
    print('no base file, quitting')
    sys.exit()

cmip6_var_max = xr.open_dataset('cmip6_output/growing_season/cmip6_%s_tasmax_grow_max_%s_%s_regrid.nc'%(crop, region, model))
cmip6_var_mean = xr.open_dataset('cmip6_output/growing_season/cmip6_%s_tasmax_grow_mean_%s_%s_regrid.nc'%(crop, region, model))

print('calc trends for %s...'%model)
tmp_var_max_trend = np.full([cmip6_var_max.lat.size, cmip6_var_max.lon.size], np.nan)
tmp_var_mean_trend = np.full([cmip6_var_mean.lat.size, cmip6_var_mean.lon.size], np.nan)

for xlat in range(tmp_var_max_trend.shape[0]):
    for ylon in range(tmp_var_max_trend.shape[1]):
        
        
        curTxx = cmip6_var_max['%s_grow_max'%var].values[:, xlat, ylon]
        nn = np.where(~np.isnan(curTxx))[0]
        if len(nn) > 10:
            X = sm.add_constant(range(len(nn)))
            mdl = sm.OLS(curTxx[nn], X).fit()
            tmp_var_max_trend[xlat, ylon] = mdl.params[1]*10

        curTmean = cmip6_var_mean['%s_grow_mean'%var].values[:, xlat, ylon]
        nn = np.where(~np.isnan(curTmean))[0]
        if len(nn) > 10:
            X = sm.add_constant(range(len(nn)))
            mdl = sm.OLS(curTmean[nn], X).fit()
            tmp_var_mean_trend[xlat, ylon] = mdl.params[1]*10

tmpDs_var_max = xr.DataArray(data   = tmp_var_max_trend, 
                  dims   = ['lat', 'lon'],
                  coords = {'lat':cmip6_var_max.lat, 'lon':cmip6_var_max.lon},
                  attrs  = {'units'     : 'C'
                    })

tmpDs_var_mean = xr.DataArray(data   = tmp_var_mean_trend, 
                  dims   = ['lat', 'lon'],
                  coords = {'lat':cmip6_var_mean.lat, 'lon':cmip6_var_mean.lon},
                  attrs  = {'units'     : 'C'
                    })

cmip6_var_max_trend = xr.Dataset()
cmip6_var_mean_trend = xr.Dataset()

cmip6_var_max_trend['%s_grow_max_trend'%var] = tmpDs_var_max.assign_coords({'model':model})
cmip6_var_mean_trend['%s_grow_mean_trend'%var] = tmpDs_var_mean.assign_coords({'model':model})

cmip6_var_max_trend.to_netcdf('cmip6_output/growing_season/cmip6_%s_%s_grow_max_trend_%s_%s_regrid.nc'%(crop, var, region, model))
cmip6_var_mean_trend.to_netcdf('cmip6_output/growing_season/cmip6_%s_%s_grow_mean_trend_%s_%s_regrid.nc'%(crop, var, region, model))