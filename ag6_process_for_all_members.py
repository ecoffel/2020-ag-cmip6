
import sys, os, glob

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


dirCmip6 = '/home/edcoffel/drive/MAX-Filer/Research/Climate-02/Data-02-edcoffel-F20/CMIP6'

model = 'access-esm1-5'
members = [x.split('/')[-2] for x in glob.glob("%s/%s/*/"%(dirCmip6, model))]

for member in members:    
    print('running %s, %s'%(model, member))
#     os.system('screen -d -m ipython ag6_extract_cmip6_grow_temp.py %s %s'%(model, member))
    os.system('screen -d -m ipython ag6_calc_yearly_grow_cmip6_bias_ensemble.py %s %s'%(model, member))
    