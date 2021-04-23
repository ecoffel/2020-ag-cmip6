import xarray as xr
import xesmf as xe
import numpy as np
import matplotlib.pyplot as plt
import scipy
import statsmodels.api as sm
import cartopy
import cartopy.crs as ccrs
import glob
import sys, os
import datetime

dirCmip6 = '/home/edcoffel/drive/MAX-Filer/Research/Climate-02/Data-02-edcoffel-F20/CMIP6'
dirERA5 = '/home/edcoffel/drive/MAX-Filer/Research/Climate-01/Data-edcoffel-F20/ERA5'
dirDeepak = '/home/edcoffel/drive/MAX-Filer/Research/Climate-01/Personal-F20/edcoffel-F20/data/projects/ag-land-climate/deepak/Maize_yield_1970_2013'

# cmip6_models = ['access-cm2', 'access-esm1-5', 'awi-cm-1-1-mr', 'bcc-csm2-mr', 'bcc-esm1', 'canesm5', 'ec-earth3', \
#                 'gfdl-cm4', 'gfdl-esm4', 'giss-e2-1-g', 'kace-1-0-g', 'fgoals-g3', 'inm-cm5-0', 'ipsl-cm6a-lr', 'miroc6', \
#                 'mpi-esm1-2-hr', 'mpi-esm1-2-lr', 'mri-esm2-0', 'noresm2-lm', 'noresm2-mm', 'sam0-unicon']

cmip6_models = ['access-cm2', 'access-esm1-5', 'awi-cm-1-1-mr', 'awi-esm-1-1-lr', 'bcc-esm1', 'canesm5', 'cmcc-esm2',
                 'ec-earth3', 'fgoals-f3-l', 'fgoals-g3', 'giss-e2-1-g', 'inm-cm4-8', 'inm-cm5-0', 'ipsl-cm6a-lr',
                 'ipsl-cm6a-lr-inca', 'kiost-esm', 'miroc6', 'mpi-esm1-2-ham', 'mpi-esm1-2-hr', 'mpi-esm1-2-lr',
                 'mri-esm2-0', 'noresm2-lm', 'noresm2-mm']

region = 'global'
var = 'tasmin'

if region == 'global':
    latRange = [-90, 90]
    lonRange = [0, 360]
elif region == 'us':
    latRange = [20, 55]
    lonRange = [220, 300]

    
def in_time_range(y):
    return (y >= 1981) & (y <= 2014)

n = 0
for i, model in enumerate(cmip6_models):
    
#     if os.path.isfile('cmip6_output/cmip6_%s_max_%s_%s.nc'%(var, region, model)):
#         continue
    
    print('opening %s...'%model)
    cmip6_temp_hist = xr.open_mfdataset('%s/%s/r1i1p1f1/historical/%s/*.nc'%(dirCmip6, model, var), concat_dim='time')
    
    print('selecting data for %s...'%model)
    cmip6_temp_hist = cmip6_temp_hist.sel(lat=slice(latRange[0], latRange[1]), \
                                             lon=slice(lonRange[0], lonRange[1]))
    cmip6_temp_hist = cmip6_temp_hist.sel(time=in_time_range(cmip6_temp_hist['time.year']))
    
    print('resampling %s...'%model)
#     cmip6_max_hist = cmip6_temp_hist.resample(time='1Y').max(dim='time')
    cmip6_monthly_max_hist = cmip6_temp_hist.resample(time='1M').max(dim='time')
#     cmip6_mean_hist = cmip6_temp_hist.resample(time='1Y').mean(dim='time')
    
    
#     cmip6_max_hist[var] -= 273.15
    cmip6_monthly_max_hist[var] -= 273.15
#     cmip6_mean_hist[var] -= 273.15
    
#     cmip6_max_ds = xr.Dataset()
    cmip6_monthly_max_ds = xr.Dataset()
#     cmip6_mean_ds = xr.Dataset()
    
    if n == 0:
#         timeVar = cmip6_max_hist.time
        timeVar_monthly = cmip6_monthly_max_hist.time
    
#     tempDs_max = xr.DataArray(data   = cmip6_max_hist[var], 
#                           dims   = ['time', 'lat', 'lon'],
#                           coords = {'time': timeVar, 'lat':cmip6_max_hist.lat, 'lon':cmip6_max_hist.lon},
#                           attrs  = {'units'     : 'C'
#                             })
#     cmip6_max_ds['%s_max'%var] = tempDs_max
    
    tempDs_monthly_max = xr.DataArray(data   = cmip6_monthly_max_hist[var], 
                          dims   = ['time', 'lat', 'lon'],
                          coords = {'time': timeVar_monthly, 'lat':cmip6_monthly_max_hist.lat, 'lon':cmip6_monthly_max_hist.lon},
                          attrs  = {'units'     : 'C'
                            })
    cmip6_monthly_max_ds['%s_monthly_max'%var] = tempDs_monthly_max
    
#     tempDs_mean = xr.DataArray(data   = cmip6_mean_hist[var], 
#                           dims   = ['time', 'lat', 'lon'],
#                           coords = {'time': timeVar, 'lat':cmip6_mean_hist.lat, 'lon':cmip6_mean_hist.lon},
#                           attrs  = {'units'     : 'C'
#                             })
#     cmip6_mean_ds['%s_mean'%var] = tempDs_mean
    
    
    print('saving netcdf...')
#     cmip6_max_ds.to_netcdf('cmip6_output/cmip6_%s_max_%s_%s9.nc'%(var, region, model))
    cmip6_monthly_max_ds.to_netcdf('cmip6_output/cmip6_%s_monthly_max_%s_%s.nc'%(var, region, model))
#     cmip6_mean_ds.to_netcdf('cmip6_output/cmip6_%s_mean_%s_%s9.nc'%(var, region, model))
    
    print()
    n += 1
# cmip6_monthly_mean_tx_ds.to_netcdf('cmip6_monthly_mean_tx_us.nc')