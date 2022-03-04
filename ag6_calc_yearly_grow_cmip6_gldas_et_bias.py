import xarray as xr
import xesmf as xe
import numpy as np
import matplotlib.pyplot as plt
import scipy
import statsmodels.api as sm
import cartopy
import cartopy.util
import cartopy.crs as ccrs
import glob
import sys, os
import pickle, gzip
import datetime

dirCmip6 = '/home/edcoffel/drive/MAX-Filer/Research/Climate-02/Data-02-edcoffel-F20/CMIP6'
dirERA5 = '/home/edcoffel/drive/MAX-Filer/Research/Climate-02/Data-02-edcoffel-F20/ERA5'
dirAgData = '/home/edcoffel/drive/MAX-Filer/Research/Climate-01/Personal-F20/edcoffel-F20/data/projects/ag-land-climate'
dirProj = '/home/edcoffel/drive/MAX-Filer/Research/Climate-01/Personal-F20/edcoffel-F20/research/2020-ag-cmip6'

with gzip.open('%s/gdd-kdd-lat-era5.dat'%dirAgData, 'rb') as f:
    era5_lat = pickle.load(f)
with gzip.open('%s/gdd-kdd-lon-era5.dat'%dirAgData, 'rb') as f:
    era5_lon = pickle.load(f)
    
crop = 'Maize'
region = 'global'
model = sys.argv[1]

tran_exist = False
canopy_evap_exist = False
soil_evap_exist = False

if os.path.isfile('cmip6_output/growing_season/cmip6_%s_grow_tran_mon_%s_%s_regrid_fixed_sh.nc'%(crop, region, model)):
    tran_exist=True
if os.path.isfile('cmip6_output/growing_season/cmip6_%s_grow_%s_mon_global_%s_regrid_fixed_sh.nc'%(crop, 'evspsblveg', model)):
    canopy_evap_exist=True
if os.path.isfile('cmip6_output/growing_season/cmip6_%s_grow_%s_mon_global_%s_regrid_fixed_sh.nc'%(crop, 'evspsblsoi', model)):
    soil_evap_exist=True

print('loading regridded et for %s'%model)
if tran_exist:
    cmip6_tran_grow = xr.open_dataset('cmip6_output/growing_season/cmip6_%s_grow_tran_mon_%s_%s_regrid_fixed_sh.nc'%(crop, region, model))

if canopy_evap_exist:
    cmip6_evap_canopy_grow = xr.open_dataset('cmip6_output/growing_season/cmip6_%s_grow_%s_mon_global_%s_regrid_fixed_sh.nc'%(crop, 'evspsblveg', model))
    
if soil_evap_exist:
    cmip6_evap_soil_grow = xr.open_dataset('cmip6_output/growing_season/cmip6_%s_grow_%s_mon_global_%s_regrid_fixed_sh.nc'%(crop, 'evspsblsoi', model))

print('loading pre-computed gldas...')
gldas_noah_tran_grow_mean_regrid = xr.open_dataset('gldas_output/gldas_%s_%s_tran_grow_mean_regrid_%s_fixed_sh.nc'%('NOAH',crop,region))
gldas_noah_soil_evap_grow_mean_regrid = xr.open_dataset('gldas_output/gldas_%s_%s_soil_evap_grow_mean_regrid_%s_fixed_sh.nc'%('NOAH', crop,region))
gldas_noah_canopy_evap_grow_mean_regrid = xr.open_dataset('gldas_output/gldas_%s_%s_canopy_evap_grow_mean_regrid_%s_fixed_sh.nc'%('NOAH',crop,region))

gldas_vic_tran_grow_mean_regrid = xr.open_dataset('gldas_output/gldas_%s_%s_tran_grow_mean_regrid_%s_fixed_sh.nc'%('VIC',crop,region))
gldas_vic_soil_evap_grow_mean_regrid = xr.open_dataset('gldas_output/gldas_%s_%s_soil_evap_grow_mean_regrid_%s_fixed_sh.nc'%('VIC', crop,region))
gldas_vic_canopy_evap_grow_mean_regrid = xr.open_dataset('gldas_output/gldas_%s_%s_canopy_evap_grow_mean_regrid_%s_fixed_sh.nc'%('VIC',crop,region))

gldas_tran_grow_mean_regrid = (gldas_noah_tran_grow_mean_regrid.evap_grow_mean.values + gldas_vic_tran_grow_mean_regrid.evap_grow_mean.values)/2
gldas_soil_evap_grow_mean_regrid = (gldas_noah_soil_evap_grow_mean_regrid.evap_grow_mean.values + gldas_vic_soil_evap_grow_mean_regrid.evap_grow_mean.values)/2
gldas_canopy_evap_grow_mean_regrid = (gldas_noah_canopy_evap_grow_mean_regrid.evap_grow_mean.values + gldas_vic_canopy_evap_grow_mean_regrid.evap_grow_mean.values)/2


yearly_tran_grow_bias = np.full([len(range(1981, 2014+1)), \
                                   cmip6_tran_grow.lat.values.shape[0], \
                                   cmip6_tran_grow.lon.values.shape[0]], np.nan)

if canopy_evap_exist:
    yearly_canopy_evap_grow_bias = np.full([len(range(1981, 2014+1)), \
                                       cmip6_evap_canopy_grow.lat.values.shape[0], \
                                       cmip6_evap_canopy_grow.lon.values.shape[0]], np.nan)

if soil_evap_exist:
    yearly_soil_evap_grow_bias = np.full([len(range(1981, 2014+1)), \
                                       cmip6_evap_soil_grow.lat.values.shape[0], \
                                       cmip6_evap_soil_grow.lon.values.shape[0]], np.nan)


print('processing %s...'%model)
for y, year in enumerate(range(1981, 2014+1)):
    print('year %d...'%year)
    for xlat in range(gldas_tran_grow_mean_regrid.shape[1]):
        for ylon in range(gldas_tran_grow_mean_regrid.shape[2]):
            
            if tran_exist:
                yearly_tran_grow_bias[y, xlat, ylon] = cmip6_tran_grow.tran_grow_mean.values[y, xlat, ylon] - \
                                                            gldas_tran_grow_mean_regrid[y, xlat, ylon]
            
            if canopy_evap_exist:
                yearly_canopy_evap_grow_bias[y, xlat, ylon] = cmip6_evap_canopy_grow.evspsblveg_grow_mean.values[y, xlat, ylon] - \
                                                            gldas_canopy_evap_grow_mean_regrid[y, xlat, ylon]
            
            if soil_evap_exist:
                yearly_soil_evap_grow_bias[y, xlat, ylon] = cmip6_evap_soil_grow.evspsblsoi_grow_mean.values[y, xlat, ylon] - \
                                                            gldas_soil_evap_grow_mean_regrid[y, xlat, ylon]
            

if tran_exist:
    with open('cmip6_output/bias/yearly-cmip6-gldas-tran-grow-bias-%s-%s-fixed-sh.dat'%(region, model), 'wb') as f:
        pickle.dump(yearly_tran_grow_bias, f)
    
if canopy_evap_exist:
    with open('cmip6_output/bias/yearly-cmip6-gldas-canopy-evap-grow-bias-%s-%s-fixed-sh.dat'%(region, model), 'wb') as f:
        pickle.dump(yearly_canopy_evap_grow_bias, f)
    
if soil_evap_exist:
    with open('cmip6_output/bias/yearly-cmip6-gldas-soil-evap-grow-bias-%s-%s-fixed-sh.dat'%(region, model), 'wb') as f:
        pickle.dump(yearly_soil_evap_grow_bias, f)