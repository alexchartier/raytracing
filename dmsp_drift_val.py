import numpy as np
import h5py
import os
import datetime as dt
import matplotlib.pyplot as plt
import nvector as nv
import pandas as pd
import pysat
import pysatMadrigal as pysatMad
import pvlib
import nc_utils
import calc_efield


def main():

    """
    Validate potential maps using DMSP cross-track drift data
    """

    # Define date range 
    starttime = dt.datetime(2019, 3, 1)
    endtime = dt.datetime(2019, 3, 2)
    timestep = dt.timedelta(minutes=2)
    lat_cutoff = 60. 

    # Define MIX potential filename format
    pot_fname_fmt = '/Users/chartat1/pymix/data/pot_sami_cond/mar19/ampere_mix_%04d-%02d-%02dT%02d-%02d-%02dZ.nc'

    # Speciy whether to download DMSP 
    dl_dmsp = False

    # Setup pysat (only needs to happen once)
    pysat.params['data_dirs'] = '/Users/chartat1/data/pysat'
    pysat.utils.registry.register(['pysatMadrigal.instruments.dmsp_ivm'])

    # Run through and process
    sats = 'f16', 'f17', 'f18'  # looks like mag not working on F15
    time = starttime
    while time < endtime:

        """ Load DMSP each day """
        if time == dt.datetime.combine(time, dt.time.min):
            dmsp = {}
            for sat in sats:
                dmsp[sat] = load_dmsp(sat, time, time + dt.timedelta(days=1), dl_dmsp=dl_dmsp)

        for sat in sats:
            if not isinstance(dmsp[sat], pd.DataFrame): 
                # no data from satellite that day
                continue

            tidx = np.abs(dmsp[sat].index - time) < timestep / 2

            if np.sum(tidx) == 0:
                # no data in time interval (e.g. sat outside of polar region)
                continue

            dmsp_sat_t = dmsp[sat][tidx]
    
            lat = dmsp_sat_t['gdlat']
            lon = dmsp_sat_t['glon']
            dirn = dmsp_sat_t['vi_dirn']

            """ Process potentials for ExB drifts at sat locs """
            pot_fname = pot_fname_fmt % (time.year, time.month, time.day, time.hour, time.minute, time.second)
            #B = 50E3 * 1E-9  # B in Tesla
            #lat = 85.
            #lon = 270.
            mix = calc_efield.calc_mix_efld(pot_fname)
            B = {
                'forward': dmsp[sat]['b_forward'], 
                'sunward': dmsp[sat]['b_perp'], 
                'down': dmsp[sat]['bd'],
            }
            B0, Be, Bn, Bu, = calc_efield.get_igrf(lat, lon, time)
            ExB_vel, ExB_dir = calc_efield.get_exb_drift(mix, B0, Bu, dirn, lat, lon)  
            breakpoint()  # get sat lat/lon
            
        time += timestep 


def load_dmsp(sat, stime, etime, dl_dmsp=False, mlat_cutoff=60, dec_rate=10):
    """ Load DMSP data using pysat. Optionally download """

    assert mlat_cutoff > 0, 'SH filtering not implemented yet'

    # Set user and password for Madrigal
    username = 'AlexChartier'
    password = 'alex.chartier@jhuapl.edu'

    # Create PySat instrument
    dmsp = pysat.Instrument(
        inst_module=pysatMad.instruments.dmsp_ivm, 
        inst_id=sat,
        clean_level='clean',
    )
    if dl_dmsp:
        dmsp.download(stime, etime, user=username, password=password)

    dmsp.load(date=stime)
    if dmsp.data.empty: 
        print('No data for %s on %s' % (sat, stime.strftime('%Y %b %d %H:%M')))
        return []

    # Process DMSP
    dmsp_vals = dmsp.data[::10] # decimate
    dmsp_vals = dmsp_vals[(dmsp_vals['mlat'] > mlat_cutoff) & np.isfinite(dmsp_vals['hor_ion_v'])]
    dmsp_times = dmsp_vals.index
    lats, lons, alts = dmsp_vals['gdlat'], dmsp_vals['glon'], dmsp_vals['gdalt']
    brng = calc_bearings(lats[:-1], lons[:-1], alts[:-1], lats[1:], lons[1:])
    dmsp_vals = dmsp_vals[:-1]
    dmsp_vals['vi_mag'] = dmsp_vals['hor_ion_v']
    dmsp_vals['vi_dirn'] = calc_vel_dirn(lats, lons, dmsp_times, brng, dmsp_vals['vi_mag'])

    return dmsp_vals


def calc_vel_dirn(lats, lons, times, brng, vels):
    """
    DMSP horizontal ion drifts provided in unintuitive format:
        HOR_ION_V: Horizontal ion vel (pos=sunward), units: m/s
    Here we figure out the direction 
    """
    vel_dirs = np.zeros(vels.shape) 
    ephem_df = pvlib.solarposition.get_solarposition(times, lats, lons)
    solaz = ephem_df['azimuth'][:-1]
    brng_solaz = zero_360(brng - solaz)
    brng_90 = zero_360(brng + 90)
    brng_270 = zero_360(brng - 90)

    # case A: sun on the right, positive velocity
    id = (brng_solaz < 180) & (vels >= 0)
    vel_dirs[id] = brng_90[id]
    
    # case B: sun on the left, positive velocity
    id = (brng_solaz >= 180) & (vels >= 0)
    vel_dirs[id] = brng_270[id]

    # case C: sun on the right, negative velocity
    id = (brng_solaz < 180) & (vels < 0)
    vel_dirs[id] = brng_90[id]
    
    # case D: sun on the left, negative velocity
    id = (brng_solaz >= 180) & (vels < 0)
    vel_dirs[id] = brng_270[id]
    
    return vel_dirs

    
def zero_360(vals):
    """ set values between 0 - 360 degrees """
    vals[vals >= 360] -= 360
    vals[vals < 0] += 360
    return vals


def calc_bearings(latas, lonas, latbs, lonbs, alts):
    """ Calculate bearings using reference ellipsoid (deg/km) """
    wgs84 = nv.FrameE(name='WGS84')

    brng_deg = np.zeros(len(latas)) * np.nan
    for ind, latb in enumerate(latbs):
        depth = -alts[ind] * 1E3
        pointB = wgs84.GeoPoint(latitude=latas[ind], longitude=lonas[ind], z=depth, degrees=True)
        lon = lonbs[ind]
        pointA = wgs84.GeoPoint(latitude=latb, longitude=lon, z=depth, degrees=True)
        p_AB_N = pointA.delta_to(pointB)  # note we want the bearing at point A
        brng_deg[ind] = p_AB_N.azimuth_deg 

    return brng_deg

    

if __name__ == '__main__':
    main()
