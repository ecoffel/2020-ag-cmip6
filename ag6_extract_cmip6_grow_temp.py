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
dirSacks = '/home/edcoffel/drive/MAX-Filer/Research/Climate-01/Personal-F20/edcoffel-F20/data/projects/ag-land-climate'


cmip6_models = ['access-cm2', 'access-esm1-5', 'awi-cm-1-1-mr', 'bcc-csm2-mr', 'bcc-esm1', 'canesm5', 'ec-earth3', \
                'gfdl-cm4', 'gfdl-esm4', 'giss-e2-1-g', 'kace-1-0-g', 'fgoals-g3', 'inm-cm5-0', 'ipsl-cm6a-lr', 'miroc6', \
                'mpi-esm1-2-hr', 'mpi-esm1-2-lr', 'mri-esm2-0', 'noresm2-lm', 'noresm2-mm', 'sam0-unicon']

region = 'global'
var = 'tasmax'
rcp = 'ssp245'
crop = 'Maize'
model = sys.argv[1]

if rcp == 'historical':
    yearRange = [1981, 2014]
else:
    yearRange = [2070, 2099]
    
if region == 'global':
    latRange = [-90, 90]
    lonRange = [0, 360]
elif region == 'us':
    latRange = [20, 55]
    lonRange = [220, 300]

sacksMaizeNc = xr.open_dataset('%s/sacks/%s.crop.calendar.fill.nc'%(dirSacks, crop))
sacksStart = sacksMaizeNc['plant'].values
sacksStart = np.roll(sacksStart, -int(sacksStart.shape[1]/2), axis=1)
sacksStart[sacksStart < 0] = np.nan
sacksEnd = sacksMaizeNc['harvest'].values
sacksEnd = np.roll(sacksEnd, -int(sacksEnd.shape[1]/2), axis=1)
sacksEnd[sacksEnd < 0] = np.nan

sacksLat = np.linspace(90, -90, 360)
sacksLon = np.linspace(0, 360, 720)
    
def in_time_range(y, y1, y2):
    return (y >= y1) & (y <= y2)


print('opening %s...'%model)
cmip6_temp_hist = xr.open_mfdataset('%s/%s/r1i1p1f1/%s/%s/*.nc'%(dirCmip6, model, rcp, var), concat_dim='time')

print('selecting data for %s...'%model)
cmip6_temp_hist = cmip6_temp_hist.sel(lat=slice(latRange[0], latRange[1]), \
                                         lon=slice(lonRange[0], lonRange[1]))
cmip6_temp_hist = cmip6_temp_hist.sel(time=in_time_range(cmip6_temp_hist['time.year'], yearRange[0], yearRange[1]))
cmip6_temp_hist.load()
cmip6_temp_hist[var] -= 273.15

# regrid sacks data to current model res
regridMesh_cur_model = xr.Dataset({'lat': (['lat'], cmip6_temp_hist.lat),
                                   'lon': (['lon'], cmip6_temp_hist.lon)})


regridder_start = xe.Regridder(xr.DataArray(data=sacksStart, dims=['lat', 'lon'], coords={'lat':sacksLat, 'lon':sacksLon}), regridMesh_cur_model, 'bilinear', reuse_weights=True)
regridder_end = xe.Regridder(xr.DataArray(data=sacksEnd, dims=['lat', 'lon'], coords={'lat':sacksLat, 'lon':sacksLon}), regridMesh_cur_model, 'bilinear', reuse_weights=True)

sacksStart_regrid = regridder_start(sacksStart)
sacksEnd_regrid = regridder_end(sacksEnd)

yearly_groups = cmip6_temp_hist.groupby('time.year').groups

yearly_grow_tmax = np.full([yearRange[1]-yearRange[0]+1, cmip6_temp_hist.lat.size, cmip6_temp_hist.lon.size], np.nan)
yearly_grow_tmean = np.full([yearRange[1]-yearRange[0]+1, cmip6_temp_hist.lat.size, cmip6_temp_hist.lon.size], np.nan)

# count up all non-nan grid cells so we can estimate percent complete
ngrid = 0
for xlat in range(cmip6_temp_hist.lat.size):
    for ylon in range(cmip6_temp_hist.lon.size):
        if ~np.isnan(sacksStart_regrid[xlat, ylon]) and ~np.isnan(sacksEnd_regrid[xlat, ylon]):
            ngrid += 1

n = 0
for xlat in range(cmip6_temp_hist.lat.size):
    
    for ylon in range(cmip6_temp_hist.lon.size):

        if ~np.isnan(sacksStart_regrid[xlat, ylon]) and ~np.isnan(sacksEnd_regrid[xlat, ylon]):
            
            if n % 100 == 0:
                print('%.2f%%'%(n/ngrid*100))

            if sacksStart_regrid[xlat, ylon] > sacksEnd_regrid[xlat, ylon]:

                # start loop on 2nd year to allow for growing season that crosses jan 1
                for y,year in enumerate(np.array(list(yearly_groups.keys()))[1:]):

                    curTmax1 = cmip6_temp_hist[var][np.array(yearly_groups[year-1])[int(sacksEnd_regrid[xlat, ylon]):], xlat, ylon]
                    curTmax2 = cmip6_temp_hist[var][np.array(yearly_groups[year])[:int(sacksStart_regrid[xlat, ylon])], xlat, ylon]

                    curTmax = np.concatenate([curTmax1, curTmax2])

                    if len(curTmax) > 0:
                        yearly_grow_tmax[y, xlat, ylon] = np.nanmax(curTmax)
                        yearly_grow_tmean[y, xlat, ylon] = np.nanmean(curTmax)
                n += 1

#                     cur_growingSeasonLen = (365-int(sacksStart_regrid[xlat, ylon])) + int(sacksEnd_regrid[xlat, ylon])

            else:

                for y,year in enumerate(np.array(list(yearly_groups.keys()))):

                    curTmax = cmip6_temp_hist[var][np.array(yearly_groups[year])[int(sacksStart_regrid[xlat, ylon]):int(sacksEnd_regrid[xlat, ylon])], xlat, ylon]
                    if len(curTmax) > 0:
                        yearly_grow_tmax[y, xlat, ylon] = np.nanmax(curTmax)
                        yearly_grow_tmean[y, xlat, ylon] = np.nanmean(curTmax)
                n += 1
#                         cur_growingSeasonLen = int(sacksEnd_regrid[xlat, ylon]) - int(sacksStart_regrid[xlat, ylon])

da_grow_tmax = xr.DataArray(data   = yearly_grow_tmax, 
                      dims   = ['time', 'lat', 'lon'],
                      coords = {'time': list(yearly_groups.keys()), 'lat':cmip6_temp_hist.lat, 'lon':cmip6_temp_hist.lon},
                      attrs  = {'units'     : 'C'
                        })
ds_grow_tmax = xr.Dataset()
ds_grow_tmax['%s_grow_max'%var] = da_grow_tmax

da_grow_tmean = xr.DataArray(data   = yearly_grow_tmean, 
                      dims   = ['time', 'lat', 'lon'],
                      coords = {'time': list(yearly_groups.keys()), 'lat':cmip6_temp_hist.lat, 'lon':cmip6_temp_hist.lon},
                      attrs  = {'units'     : 'C'
                        })
ds_grow_tmean = xr.Dataset()
ds_grow_tmean['%s_grow_mean'%var] = da_grow_tmean


print('saving netcdf...')
ds_grow_tmax.to_netcdf('cmip6_output/growing_season/cmip6_%s_grow_%s_%s_max_%s_%s.nc'%(crop, rcp, var, region, model))
ds_grow_tmean.to_netcdf('cmip6_output/growing_season/cmip6_%s_grow_%s_%s_mean_%s_%s.nc'%(crop, rcp, var, region, model))
    
    
#     print('resampling %s...'%model)
# #     cmip6_max_hist = cmip6_temp_hist.resample(time='1Y').max(dim='time')
#     cmip6_monthly_max_hist = cmip6_temp_hist.resample(time='1M').max(dim='time')
# #     cmip6_mean_hist = cmip6_temp_hist.resample(time='1Y').mean(dim='time')


# #     
#     cmip6_monthly_max_hist[var] -= 273.15
# #     cmip6_mean_hist[var] -= 273.15

# #     cmip6_max_ds = xr.Dataset()
#     cmip6_monthly_max_ds = xr.Dataset()
# #     cmip6_mean_ds = xr.Dataset()

#     if n == 0:
# #         timeVar = cmip6_max_hist.time
#         timeVar_monthly = cmip6_monthly_max_hist.time



#     tempDs_monthly_max = xr.DataArray(data   = cmip6_monthly_max_hist[var], 
#                           dims   = ['time', 'lat', 'lon'],
#                           coords = {'time': timeVar_monthly, 'lat':cmip6_monthly_max_hist.lat, 'lon':cmip6_monthly_max_hist.lon},
#                           attrs  = {'units'     : 'C'
#                             })
#     cmip6_monthly_max_ds['%s_monthly_max'%var] = tempDs_monthly_max

# #     tempDs_mean = xr.DataArray(data   = cmip6_mean_hist[var], 
# #                           dims   = ['time', 'lat', 'lon'],
# #                           coords = {'time': timeVar, 'lat':cmip6_mean_hist.lat, 'lon':cmip6_mean_hist.lon},
# #                           attrs  = {'units'     : 'C'
# #                             })
# #     cmip6_mean_ds['%s_mean'%var] = tempDs_mean


# #     cmip6_max_ds.to_netcdf('cmip6_output/cmip6_%s_max_%s_%s9.nc'%(var, region, model))
#     cmip6_monthly_max_ds.to_netcdf('cmip6_output/cmip6_%s_monthly_max_%s_%s.nc'%(var, region, model))
# #     cmip6_mean_ds.to_netcdf('cmip6_output/cmip6_%s_mean_%s_%s9.nc'%(var, region, model))
    
#     print()
#     n += 1
# cmip6_monthly_mean_tx_ds.to_netcdf('cmip6_monthly_mean_tx_us.nc')