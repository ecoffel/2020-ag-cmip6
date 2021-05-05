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
region = 'eu'

var = 'tasmin'

cmip6_var_max = xr.open_dataset('cmip6_output/cmip6_%s_max_regrid_%s_%s.nc'%(var, region, model))
cmip6_var_monthly_max = xr.open_dataset('cmip6_output/cmip6_%s_monthly_max_regrid_%s_%s.nc'%(var, region, model))
cmip6_var_mean = xr.open_dataset('cmip6_output/cmip6_%s_mean_regrid_%s_%s.nc'%(var, region, model))

print('calc trends for %s...'%model)
tmp_var_max_trend = np.full([cmip6_var_max.lat.values.shape[0], cmip6_var_max.lon.values.shape[0]], np.nan)
tmp_var_monthly_max_trend = np.full([12, cmip6_var_monthly_max.lat.values.shape[0], cmip6_var_monthly_max.lon.values.shape[0]], np.nan)
tmp_var_mean_trend = np.full([cmip6_var_mean.lat.values.shape[0], cmip6_var_mean.lon.values.shape[0]], np.nan)

monthly_groups = cmip6_var_monthly_max.groupby('time.month').groups

for xlat in range(tmp_var_max_trend.shape[0]):
    for ylon in range(tmp_var_max_trend.shape[1]):
        curTxx = cmip6_var_max['%s_max'%var].values[:, xlat, ylon]
        X = sm.add_constant(range(1981, 2015))
        mdl = sm.OLS(curTxx, X).fit()
        tmp_var_max_trend[xlat, ylon] = mdl.params[1]*10

        curTmean = cmip6_var_mean['%s_mean'%var].values[:, xlat, ylon]
        X = sm.add_constant(range(1981, 2015))
        mdl = sm.OLS(curTmean, X).fit()
        tmp_var_mean_trend[xlat, ylon] = mdl.params[1]*10

        for month in range(1, 13):
            curMonthlyTx = cmip6_var_monthly_max['%s_monthly_max'%var].values[monthly_groups[month], xlat, ylon]
            X = sm.add_constant(range(1981, 2015))
            mdl = sm.OLS(curMonthlyTx, X).fit()
            tmp_var_monthly_max_trend[month-1, xlat, ylon] = mdl.params[1]*10

tmpDs_var_max = xr.DataArray(data   = tmp_var_max_trend, 
                  dims   = ['lat', 'lon'],
                  coords = {'lat':cmip6_var_max.lat, 'lon':cmip6_var_max.lon},
                  attrs  = {'units'     : 'C'
                    })
tmpDs_var_monthly_max = xr.DataArray(data   = tmp_var_monthly_max_trend, 
                  dims   = ['month', 'lat', 'lon'],
                  coords = {'month':np.arange(1,13), 'lat':cmip6_var_monthly_max.lat, 'lon':cmip6_var_monthly_max.lon},
                  attrs  = {'units'     : 'C'
                    })
tmpDs_var_mean = xr.DataArray(data   = tmp_var_mean_trend, 
                  dims   = ['lat', 'lon'],
                  coords = {'lat':cmip6_var_mean.lat, 'lon':cmip6_var_mean.lon},
                  attrs  = {'units'     : 'C'
                    })

cmip6_var_max_trend = xr.Dataset()
cmip6_var_monthly_max_trend = xr.Dataset()
cmip6_var_mean_trend = xr.Dataset()

cmip6_var_max_trend['%s_max_trend'%var] = tmpDs_var_max.assign_coords({'model':model})
cmip6_var_monthly_max_trend['%s_monthly_max_trend'%var] = tmpDs_var_monthly_max.assign_coords({'model':model})
cmip6_var_mean_trend['%s_mean_trend'%var] = tmpDs_var_mean.assign_coords({'model':model})

cmip6_var_max_trend.to_netcdf('cmip6_output/cmip6_%s_max_trend_regrid_%s_%s.nc'%(var, region, model))
cmip6_var_monthly_max_trend.to_netcdf('cmip6_output/cmip6_%s_monthly_max_trend_regrid_%s_%s.nc'%(var, region, model))
cmip6_var_mean_trend.to_netcdf('cmip6_output/cmip6_%s_mean_trend_regrid_%s_%s.nc'%(var, region, model))