
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

for y, year in enumerate(range(1981,2019)):
    
    print('running %d'%year)
#     os.system('screen -d -m ipython ag6_extract_era5_grow_temp.py %d'%(year))
#     os.system('screen -d -m ipython ag6_extract_era5_grow_lai.py %d'%(year))
#     os.system('screen -d -m ipython ag6_extract_era5_grow_evap_from_soil.py %d'%(year))
#     os.system('screen -d -m ipython ag6_extract_era5_grow_evap_from_canopy.py %d'%(year))
#     os.system('screen -d -m ipython ag6_extract_cpc_grow_temp.py %d'%(year))
    os.system('screen -d -m ipython ag6_extract_gldas_grow_et.py %d'%(year))

    