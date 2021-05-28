import xarray as xr
import xesmf as xe
import numpy as np
import matplotlib.pyplot as plt
import scipy
import statsmodels.api as sm
import cartopy
import cartopy.crs as ccrs
import glob
import sys
import datetime

dirERA5 = '/home/edcoffel/drive/MAX-Filer/Research/Climate-02/Data-02-edcoffel-F20/ERA5'
dirSacks = '/home/edcoffel/drive/MAX-Filer/Research/Climate-01/Personal-F20/edcoffel-F20/data/projects/ag-land-climate'

crop = 'Maize'

year = int(sys.argv[1])

latRange = [-90, 90]
lonRange = [0, 360]

sacksMaizeNc = xr.open_dataset('%s/sacks/%s.crop.calendar.fill.nc'%(dirSacks, crop))
sacksStart = sacksMaizeNc['plant'].values
sacksStart = np.roll(sacksStart, -int(sacksStart.shape[1]/2), axis=1)
sacksStart[sacksStart < 0] = np.nan
sacksEnd = sacksMaizeNc['harvest'].values
sacksEnd = np.roll(sacksEnd, -int(sacksEnd.shape[1]/2), axis=1)
sacksEnd[sacksEnd < 0] = np.nan

sacksLat = np.linspace(90, -90, 360)
sacksLon = np.linspace(0, 360, 720)

regridMesh_cur_model = xr.Dataset()

n = 0

print('opening era5 %d...'%year)
slhf_era5 = xr.open_dataset('%s/daily/slhf_%d.nc'%(dirERA5, year))
slhf_era5.load()

sshf_era5 = xr.open_dataset('%s/daily/sshf_%d.nc'%(dirERA5, year))
sshf_era5.load()

print('opening era5 %d...'%(year-1))
slhf_era5_last_year = xr.open_dataset('%s/daily/slhf_%d.nc'%(dirERA5, year))
slhf_era5_last_year.load()

sshf_era5_last_year = xr.open_dataset('%s/daily/sshf_%d.nc'%(dirERA5, year))
sshf_era5_last_year.load()

slhf_era5 = slhf_era5.rename_dims(latitude='lat', longitude='lon')
slhf_era5 = slhf_era5.rename({'latitude':'lat', 'longitude':'lon'})
sshf_era5 = sshf_era5.rename_dims(latitude='lat', longitude='lon')
sshf_era5 = sshf_era5.rename({'latitude':'lat', 'longitude':'lon'})

slhf_era5_last_year = slhf_era5_last_year.rename_dims(latitude='lat', longitude='lon')
slhf_era5_last_year = slhf_era5_last_year.rename({'latitude':'lat', 'longitude':'lon'})
sshf_era5_last_year = sshf_era5_last_year.rename_dims(latitude='lat', longitude='lon')
sshf_era5_last_year = sshf_era5_last_year.rename({'latitude':'lat', 'longitude':'lon'})

# regrid sacks data to current model res
regridMesh_cur_model = xr.Dataset({'lat': (['lat'], slhf_era5.lat),
                                   'lon': (['lon'], slhf_era5.lon)})

regridder_start = xe.Regridder(xr.DataArray(data=sacksStart, dims=['lat', 'lon'], coords={'lat':sacksLat, 'lon':sacksLon}), regridMesh_cur_model, 'bilinear', reuse_weights=True)
regridder_end = xe.Regridder(xr.DataArray(data=sacksEnd, dims=['lat', 'lon'], coords={'lat':sacksLat, 'lon':sacksLon}), regridMesh_cur_model, 'bilinear', reuse_weights=True)

sacksStart_regrid = regridder_start(sacksStart)
sacksEnd_regrid = regridder_end(sacksEnd)

# count up all non-nan grid cells so we can estimate percent complete
ngrid = 0
for xlat in range(slhf_era5.lat.size):
    for ylon in range(slhf_era5.lon.size):
        if ~np.isnan(sacksStart_regrid[xlat, ylon]) and ~np.isnan(sacksEnd_regrid[xlat, ylon]):
            ngrid += 1


yearly_grow_ef = np.full([slhf_era5.lat.size, slhf_era5.lon.size], np.nan)

for xlat in range(slhf_era5.lat.size):

    for ylon in range(slhf_era5.lon.size):

        if ~np.isnan(sacksStart_regrid[xlat, ylon]) and ~np.isnan(sacksEnd_regrid[xlat, ylon]):
    
            if n % 1000 == 0:
                print('%.2f%%'%(n/ngrid*100))

            if sacksStart_regrid[xlat, ylon] > sacksEnd_regrid[xlat, ylon]:

                # start loop on 2nd year to allow for growing season that crosses jan 1
                cur_slhf_1 = slhf_era5_last_year['slhf'][int(sacksEnd_regrid[xlat, ylon]):, xlat, ylon]
                cur_slhf_2 = slhf_era5['slhf'][:int(sacksStart_regrid[xlat, ylon]), xlat, ylon]
                
                cur_sshf_1 = sshf_era5_last_year['sshf'][int(sacksEnd_regrid[xlat, ylon]):, xlat, ylon]
                cur_sshf_2 = sshf_era5['sshf'][:int(sacksStart_regrid[xlat, ylon]), xlat, ylon]

                cur_slhf = np.concatenate([cur_slhf_1, cur_slhf_2]).mean()
                cur_sshf = np.concatenate([cur_sshf_1, cur_sshf_2]).mean()

                if not np.isnan(cur_slhf):
                    ef = cur_slhf/(cur_slhf+cur_sshf)
                    if abs(ef) < 1:
                        yearly_grow_ef[xlat, ylon] = ef
                n += 1
                

            else:
                cur_slhf = slhf_era5['slhf'][int(sacksStart_regrid[xlat, ylon]):int(sacksEnd_regrid[xlat, ylon]), xlat, ylon].mean()
                cur_sshf = sshf_era5['sshf'][int(sacksStart_regrid[xlat, ylon]):int(sacksEnd_regrid[xlat, ylon]), xlat, ylon].mean()
                
                if not np.isnan(cur_slhf):
                    ef = cur_slhf/(cur_slhf+cur_sshf)
                    if abs(ef) < 1:
                        yearly_grow_ef[xlat, ylon] = ef
                n += 1

print('renaming dims...')

da_grow_ef = xr.DataArray(data   = yearly_grow_ef, 
                      dims   = ['lat', 'lon'],
                      coords = {'time':year, 'lat':temp_era5.lat, 'lon':temp_era5.lon},
                      attrs  = {'units'     : 'Fraction'
                        })
ds_grow_ef = xr.Dataset()
ds_grow_ef['ef_grow'] = da_grow_ef

print('saving netcdf...')
ds_grow_ef.to_netcdf('era5/growing_season/era5_%s_ef_grow_global_%d.nc'%(crop, year))
