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
import cartopy.crs as ccrs

import warnings
warnings.filterwarnings('ignore')

dataDirDiscovery = '/dartfs-hpc/rc/lab/C/CMIG/ecoffel/data/projects/ag-land-climate'

heatwave_pct = 90

crop = sys.argv[1]
wxData = sys.argv[2]

years = [int(sys.argv[3]), int(sys.argv[4])]

# load the sacks crop calendars

sacksMaizeNc = xr.open_dataset('%s/sacks/Maize.crop.calendar.fill.nc'%dataDirDiscovery)
sacksStart = sacksMaizeNc['plant'].values
sacksStart = np.roll(sacksStart, -int(sacksStart.shape[1]/2), axis=1)
sacksStart[sacksStart < 0] = np.nan
sacksEnd = sacksMaizeNc['harvest'].values
sacksEnd = np.roll(sacksEnd, -int(sacksEnd.shape[1]/2), axis=1)
sacksEnd[sacksEnd < 0] = np.nan

sacksLat = np.linspace(90, -90, 360)
sacksLon = np.linspace(0, 360, 720)

# days per year above heatwave_pct
heatwave_days = [] 

# for y, year in enumerate(range(years[0], years[1]+1)):
#     print('processing year %d for %s...'%(year, crop))

#     if wxData == 'cpc':
#         dsMax = xr.open_dataset('/dartfs-hpc/rc/lab/C/CMIG/CPC/tmax/tmax.%d.nc'%year, decode_cf=False)
#         dims = dsMax.dims
#         startingDate = datetime.datetime(1900, 1, 1, 0, 0, 0)
#         tDt = []

#         for curTTime in dsMax.time:
#             delta = datetime.timedelta(hours=int(curTTime.values))
#             tDt.append(startingDate + delta)
#         dsMax['time'] = tDt
        
#         dsMin = xr.open_dataset('/dartfs-hpc/rc/lab/C/CMIG/CPC/tmin/tmin.%d.nc'%year, decode_cf=False)
#         dims = dsMin.dims
#         startingDate = datetime.datetime(1900, 1, 1, 0, 0, 0)
#         tDt = []

#         for curTTime in dsMin.time:
#             delta = datetime.timedelta(hours=int(curTTime.values))
#             tDt.append(startingDate + delta)
#         dsMin['time'] = tDt
        
#         # load previous year
#         dsMaxLastYear = xr.open_dataset('/dartfs-hpc/rc/lab/C/CMIG/CPC/tmax/tmax.%d.nc'%(year-1), decode_cf=False)
        
#         dims = dsMaxLastYear.dims
#         startingDate = datetime.datetime(1900, 1, 1, 0, 0, 0)
#         tDt = []

#         for curTTime in dsMaxLastYear.time:
#             delta = datetime.timedelta(hours=int(curTTime.values))
#             tDt.append(startingDate + delta)
#         dsMaxLastYear['time'] = tDt
        
#         dsMinLastYear = xr.open_dataset('/dartfs-hpc/rc/lab/C/CMIG/CPC/tmin/tmin.%d.nc'%(year-1), decode_cf=False)
        
#         dims = dsMinLastYear.dims
#         startingDate = datetime.datetime(1900, 1, 1, 0, 0, 0)
#         tDt = []

#         for curTTime in dsMinLastYear.time:
#             delta = datetime.timedelta(hours=int(curTTime.values))
#             tDt.append(startingDate + delta)
#         dsMinLastYear['time'] = tDt
        
#     elif wxData == 'era5':
        dsMax = xr.open_dataset('/dartfs-hpc/rc/lab/C/CMIG/ERA5/daily/tasmax_%d.nc'%year, decode_cf=False)
        dims = dsMax.dims
        startingDate = datetime.datetime(year, 1, 1, 0, 0, 0)
        tDt = []

        for curTTime in dsMax.time:
            delta = datetime.timedelta(days=int(curTTime.values))
            tDt.append(startingDate + delta)
        dsMax['time'] = tDt
        dsMax['mx2t'] = dsMax['mx2t'] - 273.15
        
#         dsMin = xr.open_dataset('/dartfs-hpc/rc/lab/C/CMIG/ERA5/daily/tasmin_%d.nc'%year, decode_cf=False)
#         dims = dsMin.dims
#         startingDate = datetime.datetime(year, 1, 1, 0, 0, 0)
#         tDt = []

#         for curTTime in dsMin.time:
#             delta = datetime.timedelta(days=int(curTTime.values))
#             tDt.append(startingDate + delta)
#         dsMin['time'] = tDt
#         dsMin['mn2t'] = dsMin['mn2t'] - 273.15
        
        
        # load previous year
        dsMaxLastYear = xr.open_dataset('/dartfs-hpc/rc/lab/C/CMIG/ERA5/daily/tasmax_%d.nc'%(year-1), decode_cf=False)
        
        dims = dsMaxLastYear.dims
        startingDate = datetime.datetime(year-1, 1, 1, 0, 0, 0)
        tDt = []

        for curTTime in dsMaxLastYear.time:
            delta = datetime.timedelta(days=int(curTTime.values))
            tDt.append(startingDate + delta)
        dsMaxLastYear['time'] = tDt
        dsMaxLastYear['mx2t'] = dsMaxLastYear['mx2t'] - 273.15
        
#         dsMinLastYear = xr.open_dataset('/dartfs-hpc/rc/lab/C/CMIG/ERA5/daily/tasmin_%d.nc'%(year-1), decode_cf=False)
        
#         dims = dsMinLastYear.dims
#         startingDate = datetime.datetime(year-1, 1, 1, 0, 0, 0)
#         tDt = []

#         for curTTime in dsMinLastYear.time:
#             delta = datetime.timedelta(days=int(curTTime.values))
#             tDt.append(startingDate + delta)
#         dsMinLastYear['time'] = tDt
#         dsMinLastYear['mn2t'] = dsMinLastYear['mn2t'] - 273.15
            
    dsMax.load()
#     dsMin.load()
    dsMaxLastYear.load()
#     dsMinLastYear.load()
    
    if wxData == 'era5':
        lat = dsMax.latitude.values
        lon = dsMax.longitude.values
        tmax = dsMax.mx2t
#         tmin = dsMin.mn2t
        tmaxLastYear = dsMaxLastYear.mx2t
#         tminLastYear = dsMinLastYear.mn2t
    elif wxData == 'cpc':
        lat = dsMax.lat.values
        lon = dsMax.lon.values
        tmax = dsMax.tmax
        tmin = dsMin.tmin
        tmaxLastYear = dsMaxLastYear.tmax
        tminLastYear = dsMinLastYear.tmin
    
    if len(heatwave_days) == 0:
        heatwave_days = np.zeros([len(lat), len(lon)])
    
    for xlat in range(len(lat)):
        
        if xlat % 25 == 0:
            print('%.0f %% complete'%(xlat/len(lat)*100))
        
        for ylon in range(len(lon)):
            
            sacksNearestX = np.where((abs(sacksLat-lat[xlat]) == np.nanmin(abs(sacksLat-lat[xlat]))))[0][0]
            sacksNearestY = np.where((abs(sacksLon-lon[ylon]) == np.nanmin(abs(sacksLon-lon[ylon]))))[0][0]
            
            growingSeasonLen = 0

            if ~np.isnan(sacksStart[sacksNearestX,sacksNearestY]) and ~np.isnan(sacksEnd[sacksNearestX,sacksNearestY]):
                
                # in southern hemisphere when planting happens in fall and harvest happens in spring
                if sacksStart[sacksNearestX,sacksNearestY] > sacksEnd[sacksNearestX,sacksNearestY]:
                    curTmax = xr.concat([tmaxLastYear[int(sacksStart[sacksNearestX,sacksNearestY]):, xlat, ylon], \
                                         tmax[:int(sacksEnd[sacksNearestX,sacksNearestY]), xlat, ylon]], dim='time')
                    
                    curTmin = xr.concat([tminLastYear[int(sacksStart[sacksNearestX,sacksNearestY]):, xlat, ylon], \
                                         tmin[:int(sacksEnd[sacksNearestX,sacksNearestY]), xlat, ylon]], dim='time')
                    
                    growingSeasonLen = (365-int(sacksStart[sacksNearestX,sacksNearestY])) + int(sacksEnd[sacksNearestX,sacksNearestY])
                    
                else:
                    curTmax = tmax[int(sacksStart[sacksNearestX,sacksNearestY]):int(sacksEnd[sacksNearestX,sacksNearestY]), xlat, ylon]
                    curTmin = tmin[int(sacksStart[sacksNearestX,sacksNearestY]):int(sacksEnd[sacksNearestX,sacksNearestY]), xlat, ylon]
                    
                    growingSeasonLen = int(sacksEnd[sacksNearestX,sacksNearestY]) - int(sacksStart[sacksNearestX,sacksNearestY])
                
                # calc seasonal gdd/kdd
                curYearGdd = (curTmax.where(curTmax > t_low) + curTmin.where(curTmin > t_low))/2-t_low
                curYearKdd = curTmax.where(curTmax > t_high)-t_high
                
                # loop over weeks to get weekly kdd/gdd
                
                for w, wInd in enumerate(range(0, growingSeasonLen, 7)):
                    gddWeekly[xlat, ylon, w] = np.nansum(curYearGdd.values[wInd:wInd+7])
                    kddWeekly[xlat, ylon, w] = np.nansum(curYearKdd.values[wInd:wInd+7])
                    
                curYearGdd = curYearGdd.sum(dim='time')
                gdd[xlat, ylon] = curYearGdd.values
                
                curYearKdd = curYearKdd.sum(dim='time')
                kdd[xlat, ylon] = curYearKdd.values
    
    with gzip.open('%s/kdd-%s-%s-%d.dat'%(dataDirDiscovery, wxData, crop, year), 'wb') as f:
        pickle.dump(kdd, f)

    with gzip.open('%s/gdd-%s-%s-%d.dat'%(dataDirDiscovery, wxData, crop, year), 'wb') as f:
        pickle.dump(gdd, f)
    
#     with gzip.open('%s/kdd-weekly-%s-%s-%d.dat'%(dataDirDiscovery, wxData, crop, year), 'wb') as f:
#         pickle.dump(kddWeekly, f)

#     with gzip.open('%s/gdd-weekly-%s-%s-%d.dat'%(dataDirDiscovery, wxData, crop, year), 'wb') as f:
#         pickle.dump(gddWeekly, f)

    if not os.path.isfile('%s/gdd-kdd-lat-%s.dat'%(dataDirDiscovery, wxData)):
        with gzip.open('%s/gdd-kdd-lat-%s.dat'%(dataDirDiscovery, wxData), 'wb') as f:
            pickle.dump(lat, f)

    if not os.path.isfile('%s/gdd-kdd-lon-%s.dat'%(dataDirDiscovery, wxData)):
        with gzip.open('%s/gdd-kdd-lon-%s.dat'%(dataDirDiscovery, wxData), 'wb') as f:
            pickle.dump(lon, f)