
import sys, os
import numpy as np

# subsets = [(x,x) for x in np.arange(1980, 2018+1, 1)]

for year in range(1980, 2018+1):
    print('running %s'%str(year))
    os.system('screen -d -m ipython ag6_era5_heat_waves.py %d'%(year))