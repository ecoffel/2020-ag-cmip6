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
dirCPC = '/home/edcoffel/drive/MAX-Filer/Research/Climate-01/Data-edcoffel-F20/CPC/tmax'


region = 'global'
var = 'tmax'
crop = 'Maize'
year = int(sys.argv[1])

if region == 'global':
    latRange = [90, -90]
    lonRange = [0, 360]
elif region == 'us':
    latRange = [20, 55]
    lonRange = [220, 300]

yearRange = [1981, 2014]
    
sacksMaizeNc = xr.open_dataset('%s/sacks/%s.crop.calendar.fill.nc'%(dirSacks, crop))
sacksStart = sacksMaizeNc['plant'].values
# sacksStart = np.roll(sacksStart, -int(sacksStart.shape[1]/2), axis=1)
sacksStart[sacksStart < 0] = np.nan
sacksEnd = sacksMaizeNc['harvest'].values
# sacksEnd = np.roll(sacksEnd, -int(sacksEnd.shape[1]/2), axis=1)
sacksEnd[sacksEnd < 0] = np.nan

sacksLat = np.linspace(-90, 90, 360)
sacksLon = np.linspace(0, 360, 720)
    

print('opening cpc...')
cpc_temp_hist = xr.open_mfdataset('%s/tmax.%d.nc'%(dirCPC, year))
cpc_temp_hist_last_year = xr.open_mfdataset('%s/tmax.%d.nc'%(dirCPC, year-1))

print('selecting data...')
cpc_temp_hist = cpc_temp_hist.sel(lat=slice(latRange[0], latRange[1]), \
                                         lon=slice(lonRange[0], lonRange[1]))
cpc_temp_hist.load()

cpc_temp_hist_last_year = cpc_temp_hist_last_year.sel(lat=slice(latRange[0], latRange[1]), \
                                         lon=slice(lonRange[0], lonRange[1]))
cpc_temp_hist_last_year.load()

cpc_temp_hist = cpc_temp_hist.reindex(lat=list(reversed(cpc_temp_hist.lat)))
cpc_temp_hist_last_year = cpc_temp_hist_last_year.reindex(lat=list(reversed(cpc_temp_hist_last_year.lat)))


# regrid sacks data to current model res
regridMesh_cur_model = xr.Dataset({'lat': (['lat'], cpc_temp_hist.lat),
                                   'lon': (['lon'], cpc_temp_hist.lon)})


regridder_start = xe.Regridder(xr.DataArray(data=sacksStart, dims=['lat', 'lon'], coords={'lat':sacksLat, 'lon':sacksLon}), regridMesh_cur_model, 'bilinear', reuse_weights=True)
regridder_end = xe.Regridder(xr.DataArray(data=sacksEnd, dims=['lat', 'lon'], coords={'lat':sacksLat, 'lon':sacksLon}), regridMesh_cur_model, 'bilinear', reuse_weights=True)

sacksStart_regrid = regridder_start(sacksStart)
sacksEnd_regrid = regridder_end(sacksEnd)

# yearly_groups = cpc_temp_hist.groupby('time.year').groups

yearly_grow_temp_mean = np.full([cpc_temp_hist.lat.size, cpc_temp_hist.lon.size], np.nan)

# count up all non-nan grid cells so we can estimate percent complete
ngrid = 0
for xlat in range(cpc_temp_hist.lat.size):
    for ylon in range(cpc_temp_hist.lon.size):
        if ~np.isnan(sacksStart_regrid[xlat, ylon]) and ~np.isnan(sacksEnd_regrid[xlat, ylon]):
            ngrid += 1

           
n = 0
for xlat in range(cpc_temp_hist.lat.size):
    
    for ylon in range(cpc_temp_hist.lon.size):

        if ~np.isnan(sacksStart_regrid[xlat, ylon]) and ~np.isnan(sacksEnd_regrid[xlat, ylon]):
            
            if n % 100 == 0:
                print('%.2f%%'%(n/ngrid*100))

            if sacksStart_regrid[xlat, ylon] > sacksEnd_regrid[xlat, ylon]:

                curTemp1 = cpc_temp_hist_last_year[var][int(sacksEnd_regrid[xlat, ylon]):, xlat, ylon].values
                curTemp2 = cpc_temp_hist[var][:int(sacksStart_regrid[xlat, ylon]), xlat, ylon].values

                curTemp = np.concatenate([curTemp1, curTemp2])

                if len(curTemp) > 0:
                    yearly_grow_temp_mean[xlat, ylon] = np.nanmax(curTemp)
                n += 1

            else:
                
                curTemp = cpc_temp_hist[var][int(sacksStart_regrid[xlat, ylon]):int(sacksEnd_regrid[xlat, ylon]), xlat, ylon].values
                if len(curTemp) > 0:
                    yearly_grow_temp_mean[xlat, ylon] = np.nanmax(curTemp)
                n += 1

da_grow_temp_mean = xr.DataArray(data   = yearly_grow_temp_mean, 
                      dims   = ['lat', 'lon'],
                      coords = {'lat':cpc_temp_hist.lat, 'lon':cpc_temp_hist.lon},
                      attrs  = {'units'     : 'C'
                        })
ds_grow_temp_mean = xr.Dataset()
ds_grow_temp_mean['%s_grow_max'%var] = da_grow_temp_mean

print('saving netcdf...')
ds_grow_temp_mean.to_netcdf('cpc_output/cpc_%s_grow_max_%s_%d.nc'%(crop, region, year))
    
