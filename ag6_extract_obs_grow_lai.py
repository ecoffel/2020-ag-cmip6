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

dirSacks = '/home/edcoffel/drive/MAX-Filer/Research/Climate-01/Personal-F20/edcoffel-F20/data/projects/ag-land-climate'
dirLAI = '/home/edcoffel/drive/MAX-Filer/Research/Climate-01/Data-edcoffel-F20/LAI/ORNL-LAI-Climo'

file_var1 = 'LAI'
crop = 'Maize'


latRange = [-90, 90]
lonRange = [0, 360]


# LOADING SACKS CROP CALENDARS
sacksMaizeNc = xr.open_dataset('%s/sacks/%s.crop.calendar.fill.nc'%(dirSacks, crop))
sacksStart = sacksMaizeNc['plant'].values
sacksStart = np.roll(sacksStart, -int(sacksStart.shape[1]/2), axis=1)
sacksStart[sacksStart < 0] = np.nan
sacksEnd = sacksMaizeNc['harvest'].values
sacksEnd = np.roll(sacksEnd, -int(sacksEnd.shape[1]/2), axis=1)
sacksEnd[sacksEnd < 0] = np.nan

# THESE ARE THE LAT/LON GRIDS THAT SACKS IS ON
sacksLat = np.linspace(90, -90, 360)
sacksLon = np.linspace(0, 360, 720)

regridMesh_cur_model = xr.Dataset()

n = 0

# LOAD 1 YEAR OF ERA5 DATA
print('opening obs lai')
lai_obs = xr.open_dataset('%s/LAI_mean_monthly_1981-2015.nc4'%dirLAI)
lai_obs.load()

lai_obs = lai_obs.roll({'time':-5})

# THIS USES XESMF TO REGRID THE SACKS DATA TO ERA5 RES
# regrid sacks data to current model res
regridMesh_cur_model = xr.Dataset({'lat': (['lat'], sacksLat),
                                   'lon': (['lon'], sacksLon)})

regridder_lai_obs = xe.Regridder(xr.DataArray(data=lai_obs[file_var1], dims=['time', 'lat', 'lon'], coords={'lat':lai_obs.lat, 'lon':lai_obs.lon}), regridMesh_cur_model, 'bilinear', reuse_weights=True)

lai_obs = regridder_lai_obs(lai_obs.LAI)

# count up all non-nan grid cells so we can estimate percent complete
ngrid = 0
for xlat in range(lai_obs.lat.size):
    for ylon in range(lai_obs.lon.size):
        
        if not np.isnan(sacksStart[xlat, ylon]):
            curStart = datetime.datetime.strptime('2020%d'%(round(sacksStart[xlat, ylon])+1), '%Y%j').date().month
            sacksStart[xlat, ylon] = curStart
            
        if not np.isnan(sacksEnd[xlat, ylon]):
            curEnd = datetime.datetime.strptime('2020%d'%(round(sacksEnd[xlat, ylon])+1), '%Y%j').date().month
            sacksEnd[xlat, ylon] = curEnd
        
        if ~np.isnan(sacksStart[xlat, ylon]) and ~np.isnan(sacksEnd[xlat, ylon]):
            ngrid += 1


lai_obs_mean = np.full([lai_obs.lat.size, lai_obs.lon.size], np.nan)

# account for the shifted months - index 0 is august, so shifted to the right by 7
def cur_month(m):
    m_new = m-7
    if m_new < 0:
        m_new += 12
    return m_new
            
# THIS LOOPS OVER EVERY GRID CELL OF ERA5 AND EXTRACTS DAILY ERA5 DATA THAT FALLS WITHIN THE SACKS GROWING SEASON
# latitude loop
for xlat in range(lai_obs.lat.size):
    # longitude loop
    for ylon in range(lai_obs.lon.size):

        # if sacks calendar is defined at this grid cell
        if ~np.isnan(sacksStart[xlat, ylon]) and ~np.isnan(sacksEnd[xlat, ylon]):
    
            # just print out our progress
            if n % 1000 == 0:
                print('%.2f%%'%(n/ngrid*100))

            # there are 2 possibilities - that the planting date is before the harvest date in the current year (northern hemisphere), 
            # or that the planting date is late in the year and the harvest date is in the beginning of the next year (southern hemisphere)
            # we need to handle these two cases separately
            
            # SOUTHERN HEMISPHERE - NEED 2 YEARS OF ERA5 DATA (this year and last year)
            if sacksStart[xlat, ylon] > sacksEnd[xlat, ylon]:

                # start loop on 2nd year to allow for growing season that crosses jan 1
                cur_lai_obs1 = lai_obs[cur_month(int(sacksEnd[xlat, ylon])):, xlat, ylon]
                cur_lai_obs2 = lai_obs[:int(sacksStart[xlat, ylon]), xlat, ylon]
                
                cur_lai_obs = np.concatenate([cur_lai_obs1, cur_lai_obs2])

                if len(cur_lai_obs) > 0:
                    lai_obs_mean[xlat, ylon] = np.nanmean(cur_lai_obs)
                n += 1
                

            # NORTHERN HEMISPHERE - SIMPLER, JUST NEED 1 YEAR OF ERA5
            else:
                cur_lai_obs = lai_obs[int(sacksStart[xlat, ylon]):int(sacksEnd[xlat, ylon]), xlat, ylon]
                
                if len(cur_lai_obs) > 0:
                    lai_obs_mean[xlat, ylon] = np.nanmean(cur_lai_obs)
                    
                n += 1

print('renaming dims...')

# SAVE THE EXTRACTED TEMP DATA
da_grow_lai_mean = xr.DataArray(data   = lai_obs_mean, 
                      dims   = ['lat', 'lon'],
                      coords = {'lat':lai_obs.lat, 'lon':lai_obs.lon},
                      attrs  = {'units'     : 'LAI'
                        })
ds_grow_lai_mean = xr.Dataset()
ds_grow_lai_mean['lai_grow_mean'] = da_grow_lai_mean

print('saving netcdf...')
ds_grow_lai_mean.to_netcdf('lai_data/lai_%s_grow_mean_global.nc'%(crop))



