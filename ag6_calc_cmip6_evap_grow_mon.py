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
crop = 'Maize'
model = sys.argv[1]
var = 'evspsblsoi'

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
    
def in_time_range(y):
    return (y >= 1981) & (y <= 2014)


print('opening %s...'%model)
cmip6_evap_hist = xr.open_mfdataset('%s/%s/r1i1p1f1/historical/%s/%s_Lmon_*.nc'%(dirCmip6, model, var, var), concat_dim='time')

print('selecting data for %s...'%model)
cmip6_evap_hist = cmip6_evap_hist.sel(time=in_time_range(cmip6_evap_hist['time.year']))

cmip6_evap_hist.load()

# regrid sacks data to current model res
regridMesh_cur_model = xr.Dataset({'lat': (['lat'], cmip6_evap_hist.lat),
                                   'lon': (['lon'], cmip6_evap_hist.lon)})

regridder_start = xe.Regridder(xr.DataArray(data=sacksStart, dims=['lat', 'lon'], coords={'lat':sacksLat, 'lon':sacksLon}), regridMesh_cur_model, 'bilinear', reuse_weights=True)
regridder_end = xe.Regridder(xr.DataArray(data=sacksEnd, dims=['lat', 'lon'], coords={'lat':sacksLat, 'lon':sacksLon}), regridMesh_cur_model, 'bilinear', reuse_weights=True)

sacksStart_regrid = regridder_start(sacksStart)
sacksEnd_regrid = regridder_end(sacksEnd)

# convert sacks day to month
for xlat in range(sacksStart_regrid.shape[0]):
    for ylon in range(sacksStart_regrid.shape[1]):
        
        if not np.isnan(sacksStart_regrid[xlat, ylon]):
            curStart = datetime.datetime.strptime('2020%d'%(round(sacksStart_regrid[xlat, ylon])+1), '%Y%j').date().month
            sacksStart_regrid[xlat, ylon] = curStart-1
            
        if not np.isnan(sacksEnd_regrid[xlat, ylon]):
            curEnd = datetime.datetime.strptime('2020%d'%(round(sacksEnd_regrid[xlat, ylon])+1), '%Y%j').date().month
            sacksEnd_regrid[xlat, ylon] = curEnd-1

yearly_groups = cmip6_evap_hist.groupby('time.year').groups
yearly_grow_evap = np.full([2014-1981+1, cmip6_evap_hist.lat.size, cmip6_evap_hist.lon.size], np.nan)

# count up all non-nan grid cells so we can estimate percent complete
ngrid = 0
for xlat in range(cmip6_evap_hist.lat.size):
    for ylon in range(cmip6_evap_hist.lon.size):
        if ~np.isnan(sacksStart_regrid[xlat, ylon]) and ~np.isnan(sacksEnd_regrid[xlat, ylon]):
            ngrid += 1

n = 0
for xlat in range(cmip6_evap_hist.lat.size):
    
    for ylon in range(cmip6_evap_hist.lon.size):

        if ~np.isnan(sacksStart_regrid[xlat, ylon]) and ~np.isnan(sacksEnd_regrid[xlat, ylon]):
            
            if n % 100 == 0:
                print('%.2f%%'%(n/ngrid*100))

            if sacksStart_regrid[xlat, ylon] > sacksEnd_regrid[xlat, ylon]:

                # start loop on 2nd year to allow for growing season that crosses jan 1
                for y,year in enumerate(np.array(list(yearly_groups.keys()))[1:]):

                    cur_evap1 = cmip6_evap_hist[var][np.array(yearly_groups[year-1])[int(sacksStart_regrid[xlat, ylon]):], xlat, ylon]
                    cur_evap2 = cmip6_evap_hist[var][np.array(yearly_groups[year])[:int(sacksEnd_regrid[xlat, ylon])], xlat, ylon]

                    cur_evap = np.nanmean(np.concatenate([cur_evap1, cur_evap2]))

                    if not np.isnan(cur_evap):
                        yearly_grow_evap[y, xlat, ylon] = cur_evap*60*60*24
                        
                n += 1

            else:

                for y,year in enumerate(np.array(list(yearly_groups.keys()))):

                    cur_evap = np.nanmean(cmip6_evap_hist[var][np.array(yearly_groups[year])[int(sacksStart_regrid[xlat, ylon]):int(sacksEnd_regrid[xlat, ylon])], xlat, ylon])
                    
                    if not np.isnan(cur_evap):
                        yearly_grow_evap[y, xlat, ylon] = cur_evap*60*60*24
                n += 1

da_grow_evap = xr.DataArray(data   = yearly_grow_evap, 
                      dims   = ['time', 'lat', 'lon'],
                      coords = {'time': list(yearly_groups.keys()), 'lat':cmip6_evap_hist.lat, 'lon':cmip6_evap_hist.lon},
                      attrs  = {'units'     : 'mm'
                        })
ds_grow_evap = xr.Dataset()
ds_grow_evap['grow_evap'] = da_grow_evap

print('saving netcdf...')
ds_grow_evap.to_netcdf('cmip6_output/growing_season/cmip6_%s_grow_%s_mon_%s_%s_fixed_sh.nc'%(var, crop, region, model))
    
    