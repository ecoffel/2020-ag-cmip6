import rasterio as rio
import matplotlib.pyplot as plt 
import numpy as np
from scipy import interpolate
import statsmodels.api as sm
import scipy.stats as st
import os, sys, pickle, gzip
import datetime
import geopy.distance
import xarray as xr
import xesmf as xe
import cartopy.crs as ccrs
import glob

dirAgData = '/dartfs-hpc/rc/lab/C/CMIG/ecoffel/data/projects/ag-land-climate'
dirEra5 = '/dartfs-hpc/rc/lab/C/CMIG/ERA5'

# dirAgData = '/home/edcoffel/drive/MAX-Filer/Research/Climate-01/Personal-F20/edcoffel-F20/data/projects/ag-land-climate'
# dirEra5 = '/home/edcoffel/drive/MAX-Filer/Research/Climate-02/Data-02-edcoffel-F20/ERA5'

years = [1979, 2019]

sacksMaizeNc = xr.open_dataset('%s/sacks/Maize.crop.calendar.fill.nc'%dirAgData)
sacksStart = sacksMaizeNc['plant'].values
sacksStart = np.roll(sacksStart, -int(sacksStart.shape[1]/2), axis=1)
sacksStart[sacksStart < 0] = np.nan
sacksEnd = sacksMaizeNc['harvest'].values
sacksEnd = np.roll(sacksEnd, -int(sacksEnd.shape[1]/2), axis=1)
sacksEnd[sacksEnd < 0] = np.nan

sacksLat = np.linspace(90, -90, 360)
sacksLon = np.linspace(0, 360, 720)

era5_mx2t_quantiles = xr.open_dataset('era5_mx2t_quantiles.nc')
era5_mx2t_quantiles.load()

lat = era5_mx2t_quantiles.latitude.values
lon = era5_mx2t_quantiles.longitude.values

# regrid sacks data
regridMesh = xr.Dataset({'lat': (['lat'], lat),
                         'lon': (['lon'], lon),})

regridder_start = xe.Regridder(xr.DataArray(data=sacksStart, dims=['lat', 'lon'], coords={'lat':sacksLat, 'lon':sacksLon}), regridMesh, 'bilinear')
regridder_end = xe.Regridder(xr.DataArray(data=sacksEnd, dims=['lat', 'lon'], coords={'lat':sacksLat, 'lon':sacksLon}), regridMesh, 'bilinear')

sacksStart_regrid = regridder_start(sacksStart)
sacksEnd_regrid = regridder_end(sacksEnd)

heatwave_days = np.full([lat.size, lon.size, len(range(1979, 2018+1)), era5_mx2t_quantiles.mx2t.shape[0]], np.nan)
growing_season_len = np.full([lat.size, lon.size], np.nan)

nnLen = len(np.where(~np.isnan(np.reshape(sacksStart_regrid, [sacksStart_regrid.size, 1])))[0])

for y, year in enumerate(range(1980, 2018+1)):
    print('year %d'%year)
    dsMax = xr.open_dataset('%s/daily/tasmax_%d.nc'%(dirEra5, year))
    dsMax.load()
    dsMax['mx2t'] -= 273.15
    
    dsMaxLast = xr.open_dataset('%s/daily/tasmax_%d.nc'%(dirEra5, year-1))
    dsMaxLast.load()
    dsMaxLast['mx2t'] -= 273.15
    
    n = 0
    
    for xlat in range(len(lat)):

        for ylon in range(len(lon)):

            if ~np.isnan(sacksStart_regrid[xlat, ylon]) and ~np.isnan(sacksEnd_regrid[xlat, ylon]):

                if n % 30000 == 0:
                    print('%.0f %% complete'%(n/(nnLen)*100))

                # in southern hemisphere when planting happens in fall and harvest happens in spring
                if sacksStart_regrid[xlat, ylon] > sacksEnd_regrid[xlat, ylon]:
                    curTmax = xr.concat([dsMaxLast.mx2t[int(sacksStart_regrid[xlat, ylon]):, xlat, ylon], \
                                         dsMax.mx2t[:int(sacksEnd_regrid[xlat, ylon]), xlat, ylon]], dim='time')
                    cur_growingSeasonLen = (365-int(sacksStart_regrid[xlat, ylon])) + int(sacksEnd_regrid[xlat, ylon])
                
                else:
                    curTmax = dsMax.mx2t[int(sacksStart_regrid[xlat, ylon]):int(sacksEnd_regrid[xlat, ylon]), xlat, ylon]
                    cur_growingSeasonLen = int(sacksEnd_regrid[xlat, ylon]) - int(sacksStart_regrid[xlat, ylon])

                growing_season_len[xlat, ylon] = cur_growingSeasonLen
                for q in range(era5_mx2t_quantiles.mx2t.shape[0]):
                    heatwave_days[xlat, ylon, y, q] = np.where(curTmax.values > era5_mx2t_quantiles.mx2t[q, xlat, ylon].values)[0].size
                n += 1
    with open('era5_heat_wave_days_%d.dat'%year, 'wb') as f:
        pickle.dump(heatwave_days, f)
    
    if y == 0:
        with open('growing_season_len_maize.dat', 'wb') as f:
            pickle.dump(growing_season_len, f)
        
with open('era5_heat_wave_days.dat', 'wb') as f:
    pickle.dump(heatwave_days, f)