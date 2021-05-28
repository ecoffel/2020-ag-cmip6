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
cmip6_hfls_hist = xr.open_mfdataset('%s/%s/r1i1p1f1/historical/hfls/hfls_day_*.nc'%(dirCmip6, model), concat_dim='time')
cmip6_hfss_hist = xr.open_mfdataset('%s/%s/r1i1p1f1/historical/hfss/hfss_day_*.nc'%(dirCmip6, model), concat_dim='time')

print('selecting data for %s...'%model)
cmip6_hfls_hist = cmip6_hfls_hist.sel(time=in_time_range(cmip6_hfls_hist['time.year']))
cmip6_hfss_hist = cmip6_hfss_hist.sel(time=in_time_range(cmip6_hfss_hist['time.year']))

cmip6_hfls_hist.load()
cmip6_hfss_hist.load()


# regrid sacks data to current model res
regridMesh_cur_model = xr.Dataset({'lat': (['lat'], cmip6_hfls_hist.lat),
                                   'lon': (['lon'], cmip6_hfls_hist.lon)})

regridder_start = xe.Regridder(xr.DataArray(data=sacksStart, dims=['lat', 'lon'], coords={'lat':sacksLat, 'lon':sacksLon}), regridMesh_cur_model, 'bilinear', reuse_weights=True)
regridder_end = xe.Regridder(xr.DataArray(data=sacksEnd, dims=['lat', 'lon'], coords={'lat':sacksLat, 'lon':sacksLon}), regridMesh_cur_model, 'bilinear', reuse_weights=True)

sacksStart_regrid = regridder_start(sacksStart)
sacksEnd_regrid = regridder_end(sacksEnd)

yearly_groups = cmip6_hfls_hist.groupby('time.year').groups

yearly_grow_ef = np.full([2014-1981+1, cmip6_hfls_hist.lat.size, cmip6_hfls_hist.lon.size], np.nan)

# count up all non-nan grid cells so we can estimate percent complete
ngrid = 0
for xlat in range(cmip6_hfls_hist.lat.size):
    for ylon in range(cmip6_hfls_hist.lon.size):
        if ~np.isnan(sacksStart_regrid[xlat, ylon]) and ~np.isnan(sacksEnd_regrid[xlat, ylon]):
            ngrid += 1

n = 0
for xlat in range(cmip6_hfls_hist.lat.size):
    
    for ylon in range(cmip6_hfls_hist.lon.size):

        if ~np.isnan(sacksStart_regrid[xlat, ylon]) and ~np.isnan(sacksEnd_regrid[xlat, ylon]):
            
            if n % 100 == 0:
                print('%.2f%%'%(n/ngrid*100))

            if sacksStart_regrid[xlat, ylon] > sacksEnd_regrid[xlat, ylon]:

                # start loop on 2nd year to allow for growing season that crosses jan 1
                for y,year in enumerate(np.array(list(yearly_groups.keys()))[1:]):

                    cur_hfls1 = cmip6_hfls_hist['hfls'][np.array(yearly_groups[year-1])[int(sacksEnd_regrid[xlat, ylon]):], xlat, ylon]
                    cur_hfls2 = cmip6_hfls_hist['hfls'][np.array(yearly_groups[year])[:int(sacksStart_regrid[xlat, ylon])], xlat, ylon]
                    
                    cur_hfss1 = cmip6_hfss_hist['hfss'][np.array(yearly_groups[year-1])[int(sacksEnd_regrid[xlat, ylon]):], xlat, ylon]
                    cur_hfss2 = cmip6_hfss_hist['hfss'][np.array(yearly_groups[year])[:int(sacksStart_regrid[xlat, ylon])], xlat, ylon]

                    cur_hfls = np.nanmean(np.concatenate([cur_hfls1, cur_hfls2]))
                    cur_hfss = np.nanmean(np.concatenate([cur_hfss1, cur_hfss2]))

                    if not np.isnan(cur_hfls) and np.isnan(cur_hfss):
                        ef = cur_hfls/(cur_hfls+cur_hfss)
                        if abs(ef) <=1:
                            yearly_grow_ef[y, xlat, ylon] = ef
                        
                n += 1

            else:

                for y,year in enumerate(np.array(list(yearly_groups.keys()))):

                    cur_hfls = np.nanmean(cmip6_hfls_hist['hfls'][np.array(yearly_groups[year])[int(sacksStart_regrid[xlat, ylon]):int(sacksEnd_regrid[xlat, ylon])], xlat, ylon])
                    cur_hfss = np.nanmean(cmip6_hfss_hist['hfss'][np.array(yearly_groups[year])[int(sacksStart_regrid[xlat, ylon]):int(sacksEnd_regrid[xlat, ylon])], xlat, ylon])
                    
                    if not np.isnan(cur_hfls) and not np.isnan(cur_hfss):
                        ef = cur_hfls/(cur_hfls+cur_hfss)
                        if abs(ef) <= 1:
                            yearly_grow_ef[y, xlat, ylon] = ef
                n += 1

da_grow_ef = xr.DataArray(data   = yearly_grow_ef, 
                      dims   = ['time', 'lat', 'lon'],
                      coords = {'time': list(yearly_groups.keys()), 'lat':cmip6_hfls_hist.lat, 'lon':cmip6_hfls_hist.lon},
                      attrs  = {'units'     : 'Fraction'
                        })
ds_grow_ef = xr.Dataset()
ds_grow_ef['grow_ef'] = da_grow_ef


print('saving netcdf...')
ds_grow_ef.to_netcdf('cmip6_output/growing_season/cmip6_%s_grow_ef_%s_%s.nc'%(crop, region, model))
    
    