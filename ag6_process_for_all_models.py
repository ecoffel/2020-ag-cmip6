
import sys, os

cmip6_models_tasmax = ['access-cm2', 'access-esm1-5', 'awi-cm-1-1-mr', 'bcc-csm2-mr', 'bcc-esm1', 'canesm5', 'ec-earth3', \
                'gfdl-cm4', 'gfdl-esm4', 'giss-e2-1-g', 'kace-1-0-g', 'fgoals-g3', 'inm-cm5-0', 'ipsl-cm6a-lr', 'miroc6', \
                'mpi-esm1-2-hr', 'mpi-esm1-2-lr', 'mri-esm2-0', 'noresm2-lm', 'noresm2-mm', 'sam0-unicon']

cmip6_models_tasmin = ['access-cm2', 'access-esm1-5', 'awi-cm-1-1-mr', 'awi-esm-1-1-lr', 'bcc-esm1', 'canesm5', 'cmcc-esm2',
                 'ec-earth3', 'fgoals-f3-l', 'fgoals-g3', 'giss-e2-1-g', 'inm-cm4-8', 'inm-cm5-0', 'ipsl-cm6a-lr',
                 'ipsl-cm6a-lr-inca', 'kiost-esm', 'miroc6', 'mpi-esm1-2-ham', 'mpi-esm1-2-hr', 'mpi-esm1-2-lr',
                 'mri-esm2-0', 'noresm2-lm', 'noresm2-mm']

cmip6_models = ['access-cm2', 'access-esm1-5', 'awi-cm-1-1-mr', 'bcc-csm2-mr', 'bcc-esm1', 'canesm5', 'ec-earth3', \
                'gfdl-cm4', 'gfdl-esm4', 'giss-e2-1-g', 'kace-1-0-g', 'fgoals-g3', 'inm-cm5-0', 'ipsl-cm6a-lr', 'miroc6', \
                'mpi-esm1-2-hr', 'mpi-esm1-2-lr', 'mri-esm2-0', 'noresm2-lm', 'noresm2-mm', 'sam0-unicon']

for m, model in enumerate(cmip6_models):
    
#     if not os.path.isfile('cmip6_output/cmip6_tasmin_monthly_max_trend_regrid_global_%s.nc'%model):
    print('running %s'%model)
#     os.system('screen -d -m ipython ag6_calc_cmip6_trends.py %s'%(model))
    os.system('screen -d -m ipython ag6_calc_yearly_cmip6_bias.py %s'%(model))