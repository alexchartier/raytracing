import os
import datetime as dt
import numpy as np

"""
Script to run the raytracing analysis on SAMI3
Assumes SAMI3 output exists and WSPR HF radio link data already downloaded
Depends on the PHaRLAP package: https://www.dst.defence.gov.au/our-technologies/pharlap-provision-high-frequency-raytracing-laboratory-propagation-studies

edit the startup.m to set paths

"""

def main(time, sami_fn, wspr_fn, out_fn_fmt, hrs=np.arange(22, 24)):

    year = time.year
    month = time.month
    day = time.day 
    hours = '[%i:%i]' % (hrs[0], hrs[-1]);

    os.makedirs(os.path.dirname(out_fn_fmt), exist_ok='True')
    matlab_ex = '/Applications/MATLAB_R2022a.app/bin/matlab -batch "%s"'
    matlab_cmd = "raytrace_sami(%i, %i, %i, %s, '%s', '%s', '%s')" \
        % (year, month, day, hours, sami_fn, wspr_fn, out_fn_fmt)
    os.system(matlab_ex % matlab_cmd)


if __name__ == '__main__':

    time = dt.datetime(2014, 5, 23)
    sami_fn = time.strftime('sami3_regulargrid_elec_density_%Y%b%d.nc')
    wspr_fn = time.strftime('wsprspots-%Y-%m.csv')
    out_fn_fmt = 'links/{yyyy-mmm-dd-HHMM}.nc'
    
    main(time, sami_fn, wspr_fn, out_fn_fmt)
