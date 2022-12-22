"""
gps_pos_val.py
Calculate the GPS ionospheric positioning error based on files from MIT Haystack
"""
import numpy as np
import datetime as dt
from scipy.optimize import minimize
import scipy.interpolate
import matplotlib.pyplot as plt
import gps_utils as gu
import pandas as pd
import nc_utils



def main(sami_fn, gps_fn, out_fn, slist_fn, time):

    # Load daily
    sami = load_sami(sami_fn)
    gps = gu.load_and_preproc_mit_gps(gps_fn, slist_fn, time)

    times = np.unique(gps['time'])

    pos_errs = {
        'time': [],
        'site': [],
        'lat': [],
        'lon': [],
        'raw_pos_err': [],
        'corr_pos_err': [],
    }
    for t in times:
        time = pd.Timestamp(t)
        print(time)
        gps_t = gps.loc[gps['time']==time]
        sami_tind = sami['datetime'] == time 
        sami['int'].values = sami['dene0'][sami_tind].flatten()  # set the values for interpolation

        # Calculate model TEC corrections to the data
        corr_gps_t = proc_mod_delays(gps_t, sami)   

        # Run through selected sites and calculate position errors
        raw_pos_errs_t = calc_pos_errs(gps_t)
        corr_pos_errs_t = calc_pos_errs(corr_gps_t)
        for site, raw_pos_err in raw_pos_errs_t.items():
            corr_pos_err = corr_pos_errs_t[site]
            siterec = gps_t.loc[gps_t['gps_site'] == site].iloc[0]
            pos_errs['time'].append(time.timestamp())
            pos_errs['site'].append(site)
            pos_errs['lat'].append(siterec['gdlatr'])
            pos_errs['lon'].append(siterec['gdlonr'])
            pos_errs['raw_pos_err'].append(raw_pos_err)
            pos_errs['corr_pos_err'].append(corr_pos_err)
            
    for k, v in pos_errs.items():
        pos_errs[k] = np.array(v)

    # Store the data
    dim_defs = def_dims(pos_errs)
    var_defs = def_vars()
    nc_utils.write_nc(out_fn, var_defs, pos_errs, set_header, dim_defs)


def set_header(rootgrp, out_vars):
    rootgrp.description = 'Raw and corrected ionospheric GPS L1 3D positioning error ' +\
        'estimates, based on dual-frequency TEC observations and SAMI3 model'
    return rootgrp


def def_dims(out_vars):
    return {
        'npts': len(out_vars['time']),
    }   


def def_vars():
    return { 
        'time': {
            'units': 'Seconds since 0:00 UT 1/1/1970',
            'long_name': "Time",
            'dims': ('npts'),
            'type': 'f8',
        },
        'site': {
            'units': 'None',
            'long_name': "GPS receiver site name",
            'dims': ('npts'),
            'type': 'str',
        },
        'lat': {
            'units': 'Degrees',
            'long_name': "Geographic Latitude",
            'dims': ('npts'),
            'type': 'float',
        },
        'lon': {
            'units': 'Degrees',
            'long_name': "Geographic Longitude",
            'dims': ('npts'),
            'type': 'float',
        },
        'raw_pos_err': {
            'units': 'Metres',
            'long_name': "Ionospheric position error (3D) - uncorrected",
            'dims': ('npts'),
            'type': 'float',
        },
        'corr_pos_err': {
            'units': 'Metres',
            'long_name': "Ionospheric position error (3D) - corrected using SAMI3",
            'dims': ('npts'),
            'type': 'float',
        },
    }


def calc_pos_errs(gps, ref_rx_alt_m=0., l1tec_fac=6.13):
    """ Calculate the ionospheric positioning error based on specified GPS TEC 
    """

    sites = np.unique(gps['gps_site'])
    pos_errs = {}
    for site in sites:
        gps_s = gps.loc[gps['gps_site']==site] 
        delays = gps_s['los_tec'] / l1tec_fac
        tx_XYZ = np.stack(gps_s['tx_XYZ'].to_numpy())
        rx_XYZ = gps_s['rx_XYZ'].iloc[0]
        elm = gps_s['elm'].to_numpy()
        pos_errs[site] = get_pos_err(tx_XYZ, rx_XYZ, elm, delays)

    return pos_errs


def get_pos_err(tx, rx, elv, delays, verbose=False):
    """ Calculate positioning error. Locations are cartesian """
    assert len(elv) == len(delays), 'number of elv must match number of delays'
    range = np.zeros(delays.shape)
    for ind, delay in enumerate(delays):
        range[ind] = np.sqrt(np.sum((tx[ind, :] - rx) ** 2)) + delay

    def cost_func(rx_u):
        # Cost function based on known transmitters and ranges
        cost = 0
        for ind, rg in enumerate(range):
            cost += (np.sum((rx_u - tx[ind, :]) ** 2) - rg ** 2) ** 2 * elv[ind]
        return cost

    if np.sum(np.isnan(delays)) > 0:
        return np.nan 

    res_0 = cost_func(rx)
    res = minimize(cost_func, rx, method='Powell', options={'ftol':0.01})
    pos_err = np.sqrt(np.sum((rx - res.x) ** 2))
  
    assert pos_err < 50, 'Position error looks too big'

    if verbose:
        print('delays')
        print(*['%1.1f m, ' % d for d in delays])
        print('3D position error: %2.1f m' % pos_err)
        print('Inversion status: %s' % res.message)

    return pos_err


def proc_mod_delays(gps, sami):
    """ Calculate slant TEC from the model, and subtract off """

    # Calculate the delay
    tx = np.stack(gps['tx_XYZ'])
    rx = np.stack(gps['rx_XYZ'])
    sTEC = calc_sTEC(tx, rx, sami) / 1E16

    # Store
    corr_gps = gps.copy(deep=True) 
    corr_gps['los_tec'] -= sTEC

    return corr_gps


def calc_pos_err(out):
    """ Calculate positioning errors """
    out['pos_err'] = np.zeros(np.array(out['times']).shape) * np.nan
    for tind, time in enumerate(out['times']): 
        print(time.strftime('Calculating delay for %Y/%b/%d %H:%M'))
        try:
            out['pos_err'][tind] = get_pos_err(
                np.array(out['tx'][tind]['XYZ']), 
                np.array(out['rx'][tind]['XYZ']), 
                np.array(out['elv'][tind]), 
                out['delay'][tind],
            )
        except:
            print('Failed to process')
    return out


def calc_sTEC(tx, rx, mod):
    """
    Get the delay on each PRN
     tx: nx3 array of ECEF transmitter coords 
     rx: nx3 array of ECEF receiver coords 
    """
    npts = tx.shape[0]
    sTEC = np.zeros((npts,))
    for ind in range(npts):
        sTEC[ind] = get_slant_TEC(mod, tx[ind, :].copy(), rx[ind, :].copy())

    return sTEC


def sTEC_to_L1_delay(sTEC):
    delay = sTEC / (6.13 * 1E16)
    return delay


def get_slant_TEC(mod, tx, rx):
    """ Get the TEC from the model input file """
    stepsize = 10E3
    xi = create_raypath(tx, rx, stepsize) 
    assert np.sqrt(np.sum(mod['XYZ'].T ** 2, 1)).max() > np.sqrt(np.sum(xi ** 2, 1)).max(), 'Ray goes above model top'

    ray_ne = mod['int'](xi)
    sTEC = np.sum(ray_ne) * stepsize
    
    return sTEC


def create_raypath(tx_XYZ, rx_XYZ, stepsize=1E3, modtop=8400E3):
    xi = [] 
    dist = np.sqrt(np.sum((tx_XYZ - rx_XYZ) ** 2))
    unitv = (tx_XYZ - rx_XYZ) / dist

    assert np.sqrt(np.sum(rx_XYZ ** 2)) < modtop, 'Code assumes receiver is on the ground'
    pt = rx_XYZ.copy()
    while np.sqrt(np.sum(pt ** 2)) < modtop:
        xi.append(pt.copy())
        pt += unitv * stepsize
    xi = np.array(xi)

    return xi 


def load_sami(sami_fn):
    sami = nc_utils.ncread_vars(sami_fn)
    ALL = np.meshgrid(sami['alt'], sami['lat'], sami['lon'], indexing='ij')
    for ind, fld in enumerate(ALL):
        ALL[ind] = np.transpose(fld, (1, 2, 0))  # stupid meshgrid...

    sami['XYZ'] = gu.sphcart(ALL[0].flatten(), ALL[1].flatten(), ALL[2].flatten())
    sami['dene0'] *= 1E6
    sami['int'] = scipy.interpolate.NearestNDInterpolator(sami['XYZ'].T, sami['dene0'][0, :, :, :].flatten())
    sami['datetime'] = np.array([dt.datetime.utcfromtimestamp(t) for t in sami['time']])

    return sami


if __name__ == '__main__':

    args = 'globalDownload.py --url=http://cedar.openmadrigal.org/ --outputDir=./  --user_fullname="Alex Chartier" '+\
        '--user_email=alex.chartier@jhuapl.edu  --user_affiliation=APL --format=hdf5 --inst=8000 --kindat=3505 --startDate=03/01/2019 --endDate=03/02/2019'
   
    print("Install madrigalWeb (python3 -m pip install madrigalWeb), then ") 
    print(args)
    time = dt.datetime(2019, 3, 1, tzinfo=dt.timezone.utc)
    sami_fn_fmt = '~/data/sami3/%Y/sami3_regulargrid_elec_density_%Y%b%d.nc'
    gps_fn_fmt = '~/data/gps/mit_hdf/los_%Y%m%d.001.h5'
    out_fn_fmt = '~/data/gps_pos_errs/pos_errs_%Y%m%d.nc'
    slist_fn = '~/data/gps/sitelists/global_150.pkl'

    main(time.strftime(sami_fn_fmt), time.strftime(gps_fn_fmt), time.strftime(out_fn_fmt), slist_fn, time)

















