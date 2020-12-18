import rasterio as rio
import matplotlib.pyplot as plt 
from matplotlib.colors import Normalize
import numpy as np
import numpy.matlib
from scipy import interpolate
import statsmodels.api as sm
import statsmodels.formula.api as smf
import scipy.stats as st
import scipy
import os, sys, pickle, gzip
import datetime
import geopy.distance
import xarray as xr
import pandas as pd
import geopandas as gpd
import shapely.geometry
import cartopy
import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
import cartopy.feature as cfeature
import itertools
import random
import metpy
from metpy.plots import USCOUNTIES

import warnings
warnings.filterwarnings('ignore')

dataDirDiscovery = '/dartfs-hpc/rc/lab/C/CMIG/ecoffel/data/projects/ag-land-climate'

rebuildPointModels = False
reproject = True


crop = 'Maize'
wxData = 'era5'

yieldDataOld = False

oldStr = 'new'
if yieldDataOld:
    oldStr = 'old'

print('running analysis for %s with %s and %s deepak data'%(crop, wxData, oldStr))
    
    
def findConsec(data):
    # find longest consequtative sequence of years with yield data
    ptMax = (-1, -1)
    ptCur = (-1, -1)
    for i, val in enumerate(data):
        # start sequence
        if ~np.isnan(val) and ptCur[0] == -1:
            ptCur = (i, -1)
        #end sequence
        elif (np.isnan(val) and ptCur[0] >= 0):
            ptCur = (ptCur[0], i)
            if ptCur[1]-ptCur[0] > ptMax[1]-ptMax[0] or ptMax == (-1, -1):
                ptMax = ptCur
            ptCur = (-1, -1)
        # reached end of sequence
        elif i >= len(data)-1 and ptCur[0] >= 0:
            ptCur = (ptCur[0], i)
            if ptCur[1]-ptCur[0] > ptMax[1]-ptMax[0] or ptMax == (-1, -1):
                ptMax = ptCur
    return ptMax



if wxData == '20cr':
    tempYearRange = [1970, 2015]

    if yieldDataOld:
        yieldYearRange = [1970, 2008]
    else:
        yieldYearRange = [1970, 2013]
else:
    tempYearRange = [1981, 2018]

    if yieldDataOld:
        yieldYearRange = [1981, 2008]
    else:
        yieldYearRange = [1981, 2013]

sacksLat = np.linspace(90, -90, 360)
sacksLon = np.linspace(0, 360, 720)

# load gdd/kdd from cpc temperature data
if wxData == 'cpc':
    gdd = np.full([len(sacksLat), len(sacksLon), tempYearRange[1]-tempYearRange[0]+1], np.nan)
    kdd = np.full([len(sacksLat), len(sacksLon), tempYearRange[1]-tempYearRange[0]+1], np.nan)
elif wxData == 'era5':
    gdd = np.full([721, 1440, tempYearRange[1]-tempYearRange[0]+1], np.nan)
    kdd = np.full([721, 1440, tempYearRange[1]-tempYearRange[0]+1], np.nan)
elif wxData == '20cr':
    gdd = np.full([181, 360, tempYearRange[1]-tempYearRange[0]+1], np.nan)
    kdd = np.full([181, 360, tempYearRange[1]-tempYearRange[0]+1], np.nan)

for y, year in enumerate(range(tempYearRange[0], tempYearRange[1]+1)):
    print('loading gdd/kdd data for %d'%year)
    with gzip.open('%s/kdd-%s-%s-%d.dat'%(dataDirDiscovery, wxData, crop, year), 'rb') as f:
        curKdd = pickle.load(f)
        kdd[:, :, y] = curKdd
        
    with gzip.open('%s/gdd-%s-%s-%d.dat'%(dataDirDiscovery, wxData, crop, year), 'rb') as f:
        curGdd = pickle.load(f)
        gdd[:, :, y] = curGdd
        
with gzip.open('%s/gdd-kdd-lat-%s.dat'%(dataDirDiscovery, wxData), 'rb') as f:
    tempLat = pickle.load(f)

with gzip.open('%s/gdd-kdd-lon-%s.dat'%(dataDirDiscovery, wxData), 'rb') as f:
    tempLon = pickle.load(f)

maizeYield = []
for year in range(yieldYearRange[0],yieldYearRange[1]+1):
    
    if yieldDataOld:
        curMaizeYield = xr.open_dataset('%s/deepak/Maize_yield/Maize_areaweightedyield_%d.nc'%(dataDirDiscovery, year), decode_cf=False)
    else:
        curMaizeYield = xr.open_dataset('%s/deepak/Maize_yield_1970_2013/Maize_areaweightedyield_%d_ver12b.nc'%(dataDirDiscovery, year), decode_cf=False)
    
    curMaizeYield['time'] = [year]
    if len(maizeYield) == 0:
        maizeYield = curMaizeYield
    else:
        maizeYield = xr.concat([maizeYield, curMaizeYield], dim='time')
    
maizeYield.load()

# flip latitude axis so top is +90
if not yieldDataOld:
    latDeepak = np.flipud(maizeYield.latitude.values)
else:
    latDeepak = maizeYield.latitude.values
lonDeepak = np.roll(maizeYield.longitude.values, int(len(maizeYield.longitude)/2), axis=0)
lonDeepak[lonDeepak<0] += 360
maizeYield['Data'] = maizeYield.Data.transpose('latitude', 'longitude', 'time', 'level')
if not yieldDataOld:
    data = np.roll(np.flip(maizeYield['Data'], axis=0), int(len(lonDeepak)/2), axis=1)
else:
    data = np.roll(maizeYield['Data'], int(len(lonDeepak)/2), axis=1)

maizeYield['latitude'] = latDeepak
maizeYield['longitude'] = lonDeepak
maizeYield['Data'] = (('latitude', 'longitude', 'time'), np.squeeze(data))


if wxData == 'cpc':
    with gzip.open('%s/seasonal-precip-maize-gpcp.dat'%dataDirDiscovery, 'rb') as f:
        seasonalPrecip = pickle.load(f)
elif wxData == 'era5':
    with gzip.open('%s/seasonal-precip-maize-era5.dat'%dataDirDiscovery, 'rb') as f:
        seasonalPrecip = pickle.load(f)
elif wxData == '20cr':
    with gzip.open('%s/seasonal-precip-maize-20cr.dat'%dataDirDiscovery, 'rb') as f:
        seasonalPrecip = pickle.load(f)

        
        
    
if os.path.isfile('%s/global-point-model-data-%s-%s-%s.dat'%(dataDirDiscovery, wxData, crop, oldStr)) and not rebuildPointModels:
    with open('%s/global-point-model-data-%s-%s-%s.dat'%(dataDirDiscovery, wxData, crop, oldStr), 'rb') as f:
        print('loading saved point model data...')
        modelData = pickle.load(f)
        
        pointModels = modelData['pointModels']
        pointModels_KDD_GDD = modelData['pointModels_KDD_GDD']
        
        print('loaded %s'%('%s/global-point-model-data-%s-%s-%s.dat'%(dataDirDiscovery, wxData, crop, oldStr)))
else:
    print('building global point models...')
    
    minCropYears = 10

    pointModels = {}
    pointModels_STD = {}
    pointModels_KDD_GDD = {}
    
    for xlat in range(0, len(latDeepak)):

        if xlat % 10 == 0: print('%.0f %%'%(xlat/len(latDeepak)*100))

        pointModels[xlat] = {}
        pointModels_STD[xlat] = {}
        pointModels_KDD_GDD[xlat] = {}

        for ylon in range(0, len(lonDeepak)):

            y = maizeYield.Data.values[xlat, ylon, :]
            yNn = np.where(~np.isnan(y))[0]
#             ptMaxDeepak = findConsec(y)

            lat1 = latDeepak[xlat]
            lat2 = latDeepak[xlat] + (latDeepak[1]-latDeepak[0])

            lon1 = lonDeepak[ylon]
            lon2 = lonDeepak[ylon] + (lonDeepak[1]-lonDeepak[0])

            if len(yNn) >= minCropYears:

                indLat = [np.where(abs(tempLat-lat1) == np.nanmin(abs(tempLat-lat1)))[0][0],
                           np.where(abs(tempLat-lat2) == np.nanmin(abs(tempLat-lat2)))[0][0]]
                indLon = [np.where(abs(tempLon-lon1) == np.nanmin(abs(tempLon-lon1)))[0][0],
                           np.where(abs(tempLon-lon2) == np.nanmin(abs(tempLon-lon2)))[0][0]]

                indLatRange = np.arange(indLat[0], indLat[1]+1)
                indLonRange = np.arange(indLon[0], indLon[1]+1)

                k = np.nanmean(kdd[indLatRange, :, :], axis=0)
                k = np.nanmean(k[indLonRange, :], axis=0)

                g = np.nanmean(gdd[indLatRange, :, :], axis=0)
                g = np.nanmean(g[indLonRange, :], axis=0)
                
                indLatPr = [np.where(abs(sacksLat-lat1) == np.nanmin(abs(sacksLat-lat1)))[0][0],
                           np.where(abs(sacksLat-lat2) == np.nanmin(abs(sacksLat-lat2)))[0][0]]
                indLonPr = [np.where(abs(sacksLon-lon1) == np.nanmin(abs(sacksLon-lon1)))[0][0],
                           np.where(abs(sacksLon-lon2) == np.nanmin(abs(sacksLon-lon2)))[0][0]]

                indLatPrRange = np.arange(indLatPr[0], indLatPr[1]+1)
                indLonPrRange = np.arange(indLonPr[0], indLonPr[1]+1)

                p = np.nanmean(seasonalPrecip[indLatPrRange, :, :], axis=0)
                p = np.nanmean(p[indLonPrRange, :], axis=0)
                
                cropLen = len(y)
                k = k[:cropLen]
                g = g[:cropLen]
                p = p[:cropLen]
                
                allNn = np.where((~np.isnan(g)) & (~np.isnan(k)) & (~np.isnan(p)) & (~np.isnan(y)))[0]
                
                g = g[allNn]
                k = k[allNn]
                p = p[allNn]
                y = y[allNn]
                
                if len(np.where(np.isnan(k))[0]) == 0 and \
                    len(np.where(np.isnan(g))[0]) == 0 and \
                    len(np.where(np.isnan(p))[0]) == 0 and \
                    len(np.where(np.isnan(y))[0]) == 0:
                    
                    g = scipy.signal.detrend(g)
                    k = scipy.signal.detrend(k)
                    p = scipy.signal.detrend(p)
                    y = scipy.signal.detrend(y)
                else:
                    continue

                if len(np.where((np.isnan(k)) | (k == 0))[0]) == 0 and \
                    len(np.where((np.isnan(g)) | (g == 0))[0]) == 0 and \
                    len(np.where((np.isnan(p)) | (p == 0))[0]) == 0:

                    data = {'GDD':g, \
                            'KDD':k, \
                            'Pr':p, \
                            'Yield':y}

                    df = pd.DataFrame(data, \
                                      columns=['GDD', 'KDD', 'Pr', \
                                               'Yield'])

                    mdl = smf.ols(formula='Yield ~ KDD + GDD + Pr', data=df).fit()
                    
                    if mdl.f_pvalue <= 0.05:
                        pointModels[xlat][ylon] = mdl
                    
                    
                    
                    
                    mdl_KDD_GDD = smf.ols(formula='Yield ~ KDD + GDD', data=df).fit()
                    
                    if mdl_KDD_GDD.f_pvalue <= 0.05:
                        pointModels_KDD_GDD[xlat][ylon] = mdl_KDD_GDD
                    
                    
                    
                    
                    
                    dataStd = {'GDD':g/np.linalg.norm(g), \
                            'KDD':k/np.linalg.norm(k), \
                            'Pr':p/np.linalg.norm(p), \
                            'Yield':y/np.linalg.norm(y)}

                    dfStd = pd.DataFrame(dataStd, \
                                      columns=['GDD', 'KDD', 'Pr', \
                                               'Yield'])

                    mdlStd = smf.ols(formula='Yield ~ KDD + GDD + Pr', data=dfStd).fit()
                    
                    if mdlStd.f_pvalue <= 0.05:
                        pointModels_STD[xlat][ylon] = mdlStd
                

    with open('%s/global-point-model-data-%s-%s-%s.dat'%(dataDirDiscovery, wxData, crop, oldStr), 'wb') as f:
        modelData = {'pointModels':pointModels, \
                     'pointModels_KDD_GDD':pointModels_KDD_GDD, \
                     'pointModels_STD':pointModels_STD}
        pickle.dump(modelData, f)

        
        

if reproject:
    # project climate-related yield change with point models

    print('calculating global yield projections...')
    
    leaveOutN = 100
    
    yieldProj = np.full([len(latDeepak), len(lonDeepak), leaveOutN], np.nan)
    yieldProj_KDD_GDD = np.full([len(latDeepak), len(lonDeepak), leaveOutN], np.nan)
    globalKddChg = np.full([len(latDeepak), len(lonDeepak), leaveOutN], np.nan)
    globalGddChg = np.full([len(latDeepak), len(lonDeepak), leaveOutN], np.nan)
    globalPrChg = np.full([len(latDeepak), len(lonDeepak), leaveOutN], np.nan)

    for xlat in range(len(latDeepak)):

        if xlat not in pointModels.keys(): 
            if xlat not in pointModels_KDD_GDD.keys():
                continue
        
        if xlat % 10 == 0:
            print('%.0f %%'%(xlat/len(latDeepak)*100))

        for ylon in range(len(lonDeepak)):

            if ylon not in pointModels[xlat].keys():
                if ylon not in pointModels_KDD_GDD[xlat].keys():
                    continue
            
            lat1 = latDeepak[xlat]
            lat2 = latDeepak[xlat] + (latDeepak[1]-latDeepak[0])

            lon1 = lonDeepak[ylon]
            lon2 = lonDeepak[ylon] + (lonDeepak[1]-lonDeepak[0])

            indLat = [np.where(abs(tempLat-lat1) == np.nanmin(abs(tempLat-lat1)))[0][0],
                       np.where(abs(tempLat-lat2) == np.nanmin(abs(tempLat-lat2)))[0][0]]
            indLon = [np.where(abs(tempLon-lon1) == np.nanmin(abs(tempLon-lon1)))[0][0],
                       np.where(abs(tempLon-lon2) == np.nanmin(abs(tempLon-lon2)))[0][0]]

            indLatRange = np.arange(indLat[0], indLat[1]+1)
            indLonRange = np.arange(indLon[0], indLon[1]+1)

            indLatPr = [np.where(abs(sacksLat-lat1) == np.nanmin(abs(sacksLat-lat1)))[0][0],
                       np.where(abs(sacksLat-lat2) == np.nanmin(abs(sacksLat-lat2)))[0][0]]
            indLonPr = [np.where(abs(sacksLon-lon1) == np.nanmin(abs(sacksLon-lon1)))[0][0],
                       np.where(abs(sacksLon-lon2) == np.nanmin(abs(sacksLon-lon2)))[0][0]]

            indLatPrRange = np.arange(indLatPr[0], indLatPr[1]+1)
            indLonPrRange = np.arange(indLonPr[0], indLonPr[1]+1)

            curProjChg = []

            curKddStarts = []
            curKddEnds = []
            curGddStarts = []
            curGddEnds = []
            curPrStarts = []
            curPrEnds = []

            curKdd = np.nanmean(kdd[indLatRange, :, :], axis=0)
            curKdd = np.nanmean(curKdd[indLonRange, :], axis=0)
            
            curGdd = np.nanmean(gdd[indLatRange, :, :], axis=0)
            curGdd = np.nanmean(curGdd[indLonRange, :], axis=0)
            
            curPr = np.nanmean(seasonalPrecip[indLatPrRange, :, :], axis=0)
            curPr = np.nanmean(curPr[indLonPrRange, :], axis=0)

            gddStartLeaveOut = []
            gddEndLeaveOut = []
            
            kddStartLeaveOut = []
            kddEndLeaveOut = []
            
            prStartLeaveOut = []
            prEndLeaveOut = []
            
            for n in range(leaveOutN):
            
                inds = np.arange(0, len(curKdd))
                inds = random.sample(set(inds), len(inds)-1)
                inds.sort()
            
                X = sm.add_constant(range(len(curKdd[inds])))
                curKddMdl = sm.OLS(curKdd[inds], X).fit()
                curKddInt = curKddMdl.params[0]
                curKddTrend = curKddMdl.params[1]

                X = sm.add_constant(range(len(curGdd[inds])))
                curGddMdl = sm.OLS(curGdd[inds], X).fit()
                curGddInt = curGddMdl.params[0]
                curGddTrend = curGddMdl.params[1]

                X = sm.add_constant(range(len(curPr[inds])))
                curPrMdl = sm.OLS(curPr[inds], X).fit()
                curPrInt = curPrMdl.params[0]
                curPrTrend = curPrMdl.params[1]

                gddStartLeaveOut.append(curGddInt)
                kddStartLeaveOut.append(curKddInt)
                prStartLeaveOut.append(curPrInt)
                
                gddEndLeaveOut.append((curGddInt+curGddTrend*(2020-1979)))
                kddEndLeaveOut.append((curKddInt+curKddTrend*(2020-1979)))
                prEndLeaveOut.append((curPrInt+curPrTrend*(2020-1979)))
                
            dataStart = {'GDD':gddStartLeaveOut, \
                    'KDD':kddStartLeaveOut, \
                    'Pr':prStartLeaveOut}
            dataEnd = {'GDD':gddEndLeaveOut, \
                    'KDD':kddEndLeaveOut, \
                    'Pr':prEndLeaveOut}

            dfStart = pd.DataFrame(dataStart, columns=['GDD', 'KDD', 'Pr'])
            dfEnd = pd.DataFrame(dataEnd, columns=['GDD', 'KDD', 'Pr'])

            if xlat in pointModels.keys():
                if ylon in pointModels[xlat].keys():
                    curProjStarts = pointModels[xlat][ylon].predict(dfStart).values
                    curProjEnds = pointModels[xlat][ylon].predict(dfEnd).values

                    curProjChg = ((curProjEnds-curProjStarts)/np.nanmean(maizeYield.Data.values[xlat, ylon, :]))*100
                    curProjChg[curProjChg < -100] = np.nan
                    curProjChg[curProjChg > 100] = np.nan
                    
                    yieldProj[xlat, ylon, :] = curProjChg
            
            if xlat in pointModels_KDD_GDD.keys():
                if ylon in pointModels_KDD_GDD[xlat].keys():
            
                    curProjStarts_KDD_GDD = pointModels_KDD_GDD[xlat][ylon].predict(dfStart).values
                    curProjEnds_KDD_GDD = pointModels_KDD_GDD[xlat][ylon].predict(dfEnd).values

                    curProjChg_KDD_GDD = ((curProjEnds_KDD_GDD-curProjStarts_KDD_GDD)/np.nanmean(maizeYield.Data.values[xlat, ylon, :]))*100
                    curProjChg_KDD_GDD[curProjChg_KDD_GDD < -100] = np.nan
                    curProjChg_KDD_GDD[curProjChg_KDD_GDD > 100] = np.nan
                    
                    yieldProj_KDD_GDD[xlat, ylon, :] = curProjChg_KDD_GDD
            
            globalKddChg[xlat, ylon, :] = dfEnd['KDD'].values-dfStart['KDD'].values
            globalGddChg[xlat, ylon, :] = dfEnd['GDD'].values-dfStart['GDD'].values
            globalPrChg[xlat, ylon, :] = dfEnd['Pr'].values-dfStart['Pr'].values

            
    with open('%s/global-yield-projections-trendMethod-%s-%s-%s'%(dataDirDiscovery, crop, wxData, oldStr), 'wb') as f:
        globalYieldProj = {'yieldProj':yieldProj, \
                           'yieldProj_KDD_GDD':yieldProj_KDD_GDD, \
                           'globalKddChg':globalKddChg, \
                           'globalGddChg':globalGddChg, \
                           'globalPrChg':globalPrChg}
        pickle.dump(globalYieldProj, f)