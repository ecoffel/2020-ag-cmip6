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

dirProj = '/home/edcoffel/drive/MAX-Filer/Research/Climate-01/Personal-F20/edcoffel-F20/research/2020-ag-cmip6'

region = 'global'
model = sys.argv[1]

print('loading yearly bias for %s'%model)
with open('cmip6_output/bias/yearly-cmip6-era5-tasmax-max-bias-%s-%s.dat'%(region, model), 'rb') as f:
    yearly_tasmax_max_bias = pickle.load(f)
with open('cmip6_output/bias/yearly-cmip6-era5-tasmax-monthly-max-bias-%s-%s.dat'%(region, model), 'rb') as f:
    yearly_tasmax_monthly_max_bias = pickle.load(f)
with open('cmip6_output/bias/yearly-cmip6-era5-tasmax-mean-bias-%s-%s.dat'%(region, model), 'rb') as f:
    yearly_tasmax_mean_bias = pickle.load(f)
    
yearly_tasmax_max_bias_trend = np.full([yearly_tasmax_max_bias.shape[1], yearly_tasmax_max_bias.shape[2]], np.nan)
yearly_tasmax_max_bias_trend_sig = np.full([yearly_tasmax_max_bias.shape[1], yearly_tasmax_max_bias.shape[2]], np.nan)

yearly_tasmax_monthly_max_bias_trend = np.full([12, yearly_tasmax_max_bias.shape[1], yearly_tasmax_max_bias.shape[2]], np.nan)
yearly_tasmax_monthly_max_bias_trend_sig = np.full([12, yearly_tasmax_max_bias.shape[1], yearly_tasmax_max_bias.shape[2]], np.nan)

yearly_tasmax_mean_bias_trend = np.full([yearly_tasmax_mean_bias.shape[1], yearly_tasmax_mean_bias.shape[2]], np.nan)
yearly_tasmax_mean_bias_trend_sig = np.full([yearly_tasmax_mean_bias.shape[1], yearly_tasmax_mean_bias.shape[2]], np.nan)

print('processing trends in cmip6 bias %s...'%model)
for xlat in range(yearly_tasmax_max_bias_trend.shape[0]):
    for ylon in range(yearly_tasmax_max_bias_trend.shape[1]):
        curBias = np.squeeze(yearly_tasmax_max_bias[:, xlat, ylon])
        X = sm.add_constant(range(1981, 2015))
        mdl = sm.RLM(curBias, X).fit()
        yearly_tasmax_max_bias_trend[xlat, ylon] = mdl.params[1]*10
        yearly_tasmax_max_bias_trend_sig[xlat, ylon] = mdl.pvalues[1]

        curBias = np.squeeze(yearly_tasmax_mean_bias[:, xlat, ylon])
        X = sm.add_constant(range(1981, 2015))
        mdl = sm.RLM(curBias, X).fit()
        yearly_tasmax_mean_bias_trend[xlat, ylon] = mdl.params[1]*10
        yearly_tasmax_mean_bias_trend_sig[xlat, ylon] = mdl.pvalues[1]

        for month in range(12):
            curBias = np.squeeze(yearly_tasmax_monthly_max_bias[month, :, xlat, ylon])
            X = sm.add_constant(range(1981, 2015))
            mdl = sm.RLM(curBias, X).fit()
            yearly_tasmax_monthly_max_bias_trend[month, xlat, ylon] = mdl.params[1]*10
            yearly_tasmax_monthly_max_bias_trend_sig[month, xlat, ylon] = mdl.pvalues[1]

with open('cmip6_output/bias/yearly-cmip6-era5-tasmax-bias-trend-%s-%s.dat'%(region, model), 'wb') as f:
    pickle.dump({'yearly_tasmax_max_bias_trend': yearly_tasmax_max_bias_trend,
                 'yearly_tasmax_max_bias_trend_sig':yearly_tasmax_max_bias_trend_sig,
                 'yearly_tasmax_monthly_max_bias_trend': yearly_tasmax_monthly_max_bias_trend,
                 'yearly_tasmax_monthly_max_bias_trend_sig':yearly_tasmax_monthly_max_bias_trend_sig,
                 'yearly_tasmax_mean_bias_trend':yearly_tasmax_mean_bias_trend,
                 'yearly_tasmax_mean_bias_trend_sig':yearly_tasmax_mean_bias_trend_sig}, f)