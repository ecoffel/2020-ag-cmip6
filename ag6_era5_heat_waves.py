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

maxT = True
minT = False

crop = 'Rice'

maxVar = 'mx2t'
minVar = 'mn2t'

year = int(sys.argv[1])

# if os.path.isfile('heat_wave_days/era5_heat_wave_days_%d.dat'%year):
#     sys.exit()

sacksMaizeNc = xr.open_dataset('%s/sacks/%s.crop.calendar.fill.nc'%(dirAgData, crop))
sacksStart = sacksMaizeNc['plant'].values
sacksStart = np.roll(sacksStart, -int(sacksStart.shape[1]/2), axis=1)
sacksStart[sacksStart < 0] = np.nan
sacksEnd = sacksMaizeNc['harvest'].values
sacksEnd = np.roll(sacksEnd, -int(sacksEnd.shape[1]/2), axis=1)
sacksEnd[sacksEnd < 0] = np.nan

sacksLat = np.linspace(90, -90, 360)
sacksLon = np.linspace(0, 360, 720)

era5_max_quantiles = xr.open_dataset('era5_mx2t_quantiles.nc')
era5_max_quantiles.load()

era5_min_quantiles = xr.open_dataset('era5_mn2t_quantiles.nc')
era5_min_quantiles.load()

lat = era5_max_quantiles.latitude.values
lon = era5_max_quantiles.longitude.values

# regrid sacks data
regridMesh = xr.Dataset({'lat': (['lat'], lat),
                         'lon': (['lon'], lon),})

regridder_start = xe.Regridder(xr.DataArray(data=sacksStart, dims=['lat', 'lon'], coords={'lat':sacksLat, 'lon':sacksLon}), regridMesh, 'bilinear')
regridder_end = xe.Regridder(xr.DataArray(data=sacksEnd, dims=['lat', 'lon'], coords={'lat':sacksLat, 'lon':sacksLon}), regridMesh, 'bilinear')

sacksStart_regrid = regridder_start(sacksStart)
sacksEnd_regrid = regridder_end(sacksEnd)

heatwave_days = np.full([lat.size, lon.size, era5_max_quantiles[maxVar].shape[0]], np.nan)
coldwave_days = np.full([lat.size, lon.size, era5_max_quantiles[maxVar].shape[0]], np.nan)
growing_season_len = np.full([lat.size, lon.size], np.nan)

nnLen = len(np.where(~np.isnan(np.reshape(sacksStart_regrid, [sacksStart_regrid.size, 1])))[0])

print('year %d'%year)
if maxT:
    if maxVar == 'mx2t':
        dsMax = xr.open_dataset('%s/daily/tasmax_%d.nc'%(dirEra5, year))
        dsMax.load()

        dsMaxLast = xr.open_dataset('%s/daily/tasmax_%d.nc'%(dirEra5, year-1))
        dsMaxLast.load()
        
        dsMax['mx2t'] -= 273.15
        dsMaxLast['mx2t'] -= 273.15
    elif maxVar == 'tw':
        dsMax = xr.open_dataset('%s/daily/tw_max_%d.nc'%(dirEra5, year))
        dsMax.load()

        dsMaxLast = xr.open_dataset('%s/daily/tw_max_%d.nc'%(dirEra5, year-1))
        dsMaxLast.load()

if minT:
    
    if minVar == 'mn2t':
        dsMin = xr.open_dataset('%s/daily/tasmin_%d.nc'%(dirEra5, year))
        dsMin.load()

        dsMinLast = xr.open_dataset('%s/daily/tasmin_%d.nc'%(dirEra5, year-1))
        dsMinLast.load()
        
        dsMin['mn2t'] -= 273.15
        dsMinLast['mn2t'] -= 273.15
    elif minVar == 'tw':
        dsMin = xr.open_dataset('%s/daily/tw_min_%d.nc'%(dirEra5, year))
        dsMin.load()
    
        dsMinLast = xr.open_dataset('%s/daily/tw_min_%d.nc'%(dirEra5, year-1))
        dsMinLast.load()
    

n = 0

for xlat in range(len(lat)):

    for ylon in range(len(lon)):

        if ~np.isnan(sacksStart_regrid[xlat, ylon]) and ~np.isnan(sacksEnd_regrid[xlat, ylon]):

            if n % 30000 == 0:
                print('%.0f %% complete'%(n/(nnLen)*100))

            # in southern hemisphere when planting happens in fall and harvest happens in spring
            if sacksStart_regrid[xlat, ylon] > sacksEnd_regrid[xlat, ylon]:
                if maxT:
                    curTmax = xr.concat([dsMaxLast[maxVar][int(sacksStart_regrid[xlat, ylon]):, xlat, ylon], \
                                         dsMax[maxVar][:int(sacksEnd_regrid[xlat, ylon]), xlat, ylon]], dim='time')
                if minT:
                    curTmin = xr.concat([dsMinLast[minVar][int(sacksStart_regrid[xlat, ylon]):, xlat, ylon], \
                                         dsMin[minVar][:int(sacksEnd_regrid[xlat, ylon]), xlat, ylon]], dim='time')
                cur_growingSeasonLen = (365-int(sacksStart_regrid[xlat, ylon])) + int(sacksEnd_regrid[xlat, ylon])

            else:
                if maxT:
                    curTmax = dsMax[maxVar][int(sacksStart_regrid[xlat, ylon]):int(sacksEnd_regrid[xlat, ylon]), xlat, ylon]
                if minT:
                    curTmin = dsMin[minVar][int(sacksStart_regrid[xlat, ylon]):int(sacksEnd_regrid[xlat, ylon]), xlat, ylon]
                cur_growingSeasonLen = int(sacksEnd_regrid[xlat, ylon]) - int(sacksStart_regrid[xlat, ylon])

            growing_season_len[xlat, ylon] = cur_growingSeasonLen
            if maxT:
                for q in range(era5_max_quantiles[maxVar].shape[0]):
                    if q < 3:
                        heatwave_days[xlat, ylon, q] = np.where(curTmax.values < era5_max_quantiles[maxVar][q, xlat, ylon].values)[0].size
                    else:
                        heatwave_days[xlat, ylon, q] = np.where(curTmax.values > era5_max_quantiles[maxVar][q, xlat, ylon].values)[0].size
            
            if minT:
                for q in range(era5_min_quantiles[minVar].shape[0]):
                    if q < 3:
                        coldwave_days[xlat, ylon, q] = np.where(curTmin.values < era5_min_quantiles[minVar][q, xlat, ylon].values)[0].size
                    else:
                        coldwave_days[xlat, ylon, q] = np.where(curTmin.values > era5_min_quantiles[minVar][q, xlat, ylon].values)[0].size
            n += 1
if maxT:
    with open('heat_wave_days/era5_heat_wave_days_%s_%d.dat'%(crop, year), 'wb') as f:
        pickle.dump(heatwave_days, f)
if minT:
    with open('heat_wave_days/era5_cold_wave_days_%s_%d.dat'%(crop, year), 'wb') as f:
        pickle.dump(coldwave_days, f)

if year == 1980:
    with open('heat_wave_days/growing_season_len_maize_%s.dat'%crop, 'wb') as f:
        pickle.dump(growing_season_len, f)
        
# with open('era5_heat_wave_days.dat', 'wb') as f:
#     pickle.dump(heatwave_days, f)