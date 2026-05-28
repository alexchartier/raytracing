import numpy as np
import datetime as dt
import nc_utils
import os

time = dt.datetime(2019, 3, 1)
endtime = dt.datetime(2019, 3, 31)
timestep = dt.timedelta(minutes=2)
pot_fname_fmt = 'ampere_mix_%04d-%02d-%02dT%02d-%02d-%02dZ.nc'
#dirname_fmt = '/disks/d0390/project/nimo/steelrj1/swo2r/pipeline/data/sami-ampere/%Y%m%d/pymixout/'
dirname_fmt = '/Users/chartat1/pymix/data/pot_sami_cond/mar19/'

while time < endtime:
    dirname = time.strftime(dirname_fmt)
    pot_fname = pot_fname_fmt % (time.year, time.month, time.day, time.hour, time.minute, time.second)
    try:
        pot = nc_utils.ncread_vars(os.path.join(dirname, pot_fname))
        crosscap = pot['Potential'].max() - pot['Potential'].min() 
        print('%s ::: %2.2f' % (time.strftime('%Y%b%d %H:%M'), crosscap))
    except:
        #print('%s not found' % pot_fname)
        None

    time += timestep

    

