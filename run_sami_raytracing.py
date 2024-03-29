import os
import datetime as dt
import numpy as np

"""
Script to run the raytracing analysis on SAMI3
Assumes SAMI3 output exists and WSPR HF radio link data already downloaded
Depends on the PHaRLAP package: https://www.dst.defence.gov.au/our-technologies/pharlap-provision-high-frequency-raytracing-laboratory-propagation-studies

edit the startup.m to set paths

"""

def main(stime, etime, sami_fn_fmt, wspr_fn_fmt, out_fn_fmt, wspr_out_fn_fmt, hrs=np.arange(0, 24)):

    time = stime

    while time <= etime:
        print(time.strftime('Processing %Y %b %d'))
        sami_fn = os.path.expanduser(time.strftime(sami_fn_fmt))
        wspr_fn = os.path.expanduser(time.strftime(wspr_fn_fmt))

        year = time.year
        month = time.month
        day = time.day 
        hours = '[%i:%i]' % (hrs[0], hrs[-1])

        os.makedirs(os.path.dirname(out_fn_fmt), exist_ok='True')
        os.makedirs(os.path.dirname(wspr_out_fn_fmt), exist_ok='True')
    #    matlab_ex = '/Applications/MATLAB_R2022a.app/bin/matlab -batch "%s"'
        matlab_ex = '/Applications/MATLAB_R2022a.app/bin/matlab -nosplash -nodisplay -nodesktop -batch "%s"'
        matlab_cmd = "raytrace_sami(%i, %i, %i, %s, '%s', '%s', '%s', '%s')" \
            % (year, month, day, hours, sami_fn, wspr_fn, out_fn_fmt, wspr_out_fn_fmt)
        print(matlab_cmd)
        os.system(matlab_ex % matlab_cmd)

        time += dt.timedelta(days=1)


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

    stime = dt.datetime(2019, 3, 30)
    etime = dt.datetime(2019, 3, 30)
    sami_fn_fmt = '~/data/sami3/%Y/sami3_regulargrid_elec_density_%Y%b%d.nc'
    wspr_fn_fmt = 'data/wspr/wsprspots-%Y-%m.csv'
    wspr_out_fn_fmt = 'data/wspr_out/wspr_{yyyy-mmm-dd-HHMM}.mat'
    out_fn_fmt = 'links/{yyyy-mmm-dd-HHMM}.nc'
    
    main(stime, etime, sami_fn_fmt, wspr_fn_fmt, out_fn_fmt, wspr_out_fn_fmt)




