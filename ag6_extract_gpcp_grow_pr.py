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
dirGPCP = '/home/edcoffel/drive/MAX-Filer/Research/Climate-01/Data-edcoffel-F20/GPCP'


region = 'global'
var = 'precip'
crop = 'Maize'

if region == 'global':
    latRange = [-90, 90]
    lonRange = [0, 360]
elif region == 'us':
    latRange = [20, 55]
    lonRange = [220, 300]

yearRange = [1981, 2014]
    
sacksMaizeNc = xr.open_dataset('%s/sacks/%s.crop.calendar.fill.nc'%(dirSacks, crop))
sacksStart = sacksMaizeNc['plant'].values
sacksStart = np.roll(sacksStart, -int(sacksStart.shape[1]/2), axis=1)
sacksStart[sacksStart < 0] = np.nan
sacksEnd = sacksMaizeNc['harvest'].values
sacksEnd = np.roll(sacksEnd, -int(sacksEnd.shape[1]/2), axis=1)
sacksEnd[sacksEnd < 0] = np.nan

sacksLat = np.linspace(90, -90, 360)
sacksLon = np.linspace(0, 360, 720)
    
def in_time_range(y):
    return (y >= yearRange[0]) & (y <= yearRange[1])


print('opening gpcp...')
gpcp_pr_hist = xr.open_dataset('%s/precip.mon.mean.nc'%dirGPCP)

print('selecting data...')
gpcp_pr_hist = gpcp_pr_hist.sel(lat=slice(latRange[0], latRange[1]), \
                                         lon=slice(lonRange[0], lonRange[1]))
gpcp_pr_hist = gpcp_pr_hist.sel(time=in_time_range(gpcp_pr_hist['time.year']))
gpcp_pr_hist.load()

# regrid sacks data to current model res
regridMesh_cur_model = xr.Dataset({'lat': (['lat'], gpcp_pr_hist.lat),
                                   'lon': (['lon'], gpcp_pr_hist.lon)})


regridder_start = xe.Regridder(xr.DataArray(data=sacksStart, dims=['lat', 'lon'], coords={'lat':sacksLat, 'lon':sacksLon}), regridMesh_cur_model, 'bilinear', reuse_weights=True)
regridder_end = xe.Regridder(xr.DataArray(data=sacksEnd, dims=['lat', 'lon'], coords={'lat':sacksLat, 'lon':sacksLon}), regridMesh_cur_model, 'bilinear', reuse_weights=True)

sacksStart_regrid = regridder_start(sacksStart)
sacksEnd_regrid = regridder_end(sacksEnd)

yearly_groups = gpcp_pr_hist.groupby('time.year').groups

yearly_grow_pr_mean = np.full([yearRange[1]-yearRange[0]+1, gpcp_pr_hist.lat.size, gpcp_pr_hist.lon.size], np.nan)

# count up all non-nan grid cells so we can estimate percent complete
ngrid = 0
for xlat in range(gpcp_pr_hist.lat.size):
    for ylon in range(gpcp_pr_hist.lon.size):
        if ~np.isnan(sacksStart_regrid[xlat, ylon]) and ~np.isnan(sacksEnd_regrid[xlat, ylon]):
            ngrid += 1

            
# convert sacks day to month
for xlat in range(sacksStart_regrid.shape[0]):
    for ylon in range(sacksStart_regrid.shape[1]):
        
        if not np.isnan(sacksStart_regrid[xlat, ylon]):
            curStart = datetime.datetime.strptime('2020%d'%(round(sacksStart_regrid[xlat, ylon])+1), '%Y%j').date().month
            sacksStart_regrid[xlat, ylon] = curStart
            
        if not np.isnan(sacksEnd_regrid[xlat, ylon]):
            curEnd = datetime.datetime.strptime('2020%d'%(round(sacksEnd_regrid[xlat, ylon])+1), '%Y%j').date().month
            sacksEnd_regrid[xlat, ylon] = curEnd
            
n = 0
for xlat in range(gpcp_pr_hist.lat.size):
    
    for ylon in range(gpcp_pr_hist.lon.size):

        if ~np.isnan(sacksStart_regrid[xlat, ylon]) and ~np.isnan(sacksEnd_regrid[xlat, ylon]):
            
            if n % 100 == 0:
                print('%.2f%%'%(n/ngrid*100))

            if sacksStart_regrid[xlat, ylon] > sacksEnd_regrid[xlat, ylon]:

                # start loop on 2nd year to allow for growing season that crosses jan 1
                for y,year in enumerate(np.array(list(yearly_groups.keys()))[1:]):

                    curPr1 = gpcp_pr_hist[var][np.array(yearly_groups[year-1])[int(sacksEnd_regrid[xlat, ylon]):], xlat, ylon].values
                    curPr2 = gpcp_pr_hist[var][np.array(yearly_groups[year])[:int(sacksStart_regrid[xlat, ylon])], xlat, ylon].values

                    curPr = np.concatenate([curPr1, curPr2])

                    if len(curPr) > 0:
                        yearly_grow_pr_mean[y, xlat, ylon] = np.nanmean(curPr)
                n += 1

            else:
                
                for y,year in enumerate(np.array(list(yearly_groups.keys()))):
                
                    curPr = gpcp_pr_hist[var][np.array(yearly_groups[year])[int(sacksStart_regrid[xlat, ylon]):int(sacksEnd_regrid[xlat, ylon])], xlat, ylon].values
                    if len(curPr) > 0:
                        yearly_grow_pr_mean[y, xlat, ylon] = np.nanmean(curPr)
                n += 1

da_grow_pr_mean = xr.DataArray(data   = yearly_grow_pr_mean, 
                      dims   = ['time', 'lat', 'lon'],
                      coords = {'time': list(yearly_groups.keys()), 'lat':gpcp_pr_hist.lat, 'lon':gpcp_pr_hist.lon},
                      attrs  = {'units'     : 'mm'
                        })
ds_grow_pr_mean = xr.Dataset()
ds_grow_pr_mean['%s_grow_mean'%var] = da_grow_pr_mean

print('saving netcdf...')
ds_grow_pr_mean.to_netcdf('gpcp_output/gpcp_%s_grow_mean_%s.nc'%(crop, region))
    
