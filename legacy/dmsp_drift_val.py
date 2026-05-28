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
from ppigrf import igrf
import scipy.interpolate


def main(
    stime=dt.datetime(2019, 3, 2, 10, 0),
    etime=dt.datetime(2019, 3, 3),
    timestep=dt.timedelta(minutes=2),
    mlat_cutoff=60.,
    pot_fname_fmt='/Users/chartat1/pymix/data/pot_sami_cond/mar19/ampere_mix_%04d-%02d-%02dT%02d-%02d-%02dZ.nc',
    plot_fname_fmt='./plots/ampere_mix_%04d-%02d-%02dT%02d-%02d-%02dZ.png',
    dl_dmsp=False,
    sats = ['f16', 'f17', 'f18'],  # looks like mag not working on F15
    dmsp_dec_rate=10,
    pysat_dir='/Users/chartat1/data/pysat',
):

    """
    Validate potential maps using DMSP cross-track drift data
    """

    """ Setup pysat (only needs to happen once) """
    pysat.params['data_dirs'] = pysat_dir
    pysat.utils.registry.register(['pysatMadrigal.instruments.dmsp_ivm'])

    """ Call IGRF for the AMPERE grid (only need this once per run too) """
    pot_fname = pot_fname_fmt % (stime.year, stime.month, stime.day, stime.hour, stime.minute, stime.second)
    pot = nc_utils.ncread_vars(pot_fname)
    mix = calc_mix_efld(pot)  
    B0, Be, Bn, Bu, = load_igrf(mix['glat'], mix['glon'], stime)

    """ Get the North Magnetic Pole location in geographic coordinates """
    np_idx = pot['MLAT (AACGM)'] == 90
    np_latlon = [pot['Geographic Latitude'][np_idx][0], pot['Geographic Longitude'][np_idx][0]]

    """ Run through and process """
    time = stime
    while time < etime:

        """ Load DMSP each day """
        if (time == dt.datetime.combine(time, dt.time.min)) or (time == stime):
            dmsp = {}
            for sat in sats:
                dmsp_data = load_dmsp(
                    sat, time, time + dt.timedelta(days=1), 
                    dl_dmsp=dl_dmsp, 
                )
                dmsp[sat] = proc_dmsp(
                    dmsp_data, 
                    dec_rate=dmsp_dec_rate, mlat_cutoff=mlat_cutoff, np_latlon=np_latlon,
                )
            

        for sat in sats:
            if not isinstance(dmsp[sat], pd.DataFrame): 
                # no data from satellite that day
                continue

            #tidx = np.abs(dmsp[sat].index - time) < timestep / 2
            tidx = np.abs(dmsp[sat].index - time) < dt.timedelta(minutes=10) / 2

            if np.sum(tidx) == 0:
                # no data in time interval (e.g. sat outside of polar region)
                continue

            """ Process potentials for ExB drifts at AMPERE grid locs """
            pot_fname = pot_fname_fmt % (time.year, time.month, time.day, time.hour, time.minute, time.second)
            pot = nc_utils.ncread_vars(pot_fname)
            mix = calc_mix_efld(pot)  
            ExB_vel, ExB_dir, vi_E, vi_N = calc_exb_drift(mix, B0, Bu, mix['glat'], mix['glon'])  


            """ Interpolate model ExB drifts to DMSP locations """
            

        
            """  plot """
            plot_fname = plot_fname_fmt % (time.year, time.month, time.day, time.hour, time.minute, time.second)
            fig, ax = plt.subplots(1, 1,subplot_kw={'projection': 'polar'}) 
    
            rad, theta = get_rad_theta(
                np.unique(pot['MLAT (AACGM)']), np.unique(pot['MLON (AACGM)']),
            )

            """ contour data over the map. """
            latlim = 55.
            im = ax.contourf(theta, rad, pot['Potential'][::-1, :-1])
            yticks = 90 - np.rad2deg(ax.get_yticks())
            ax.set_yticklabels(['%1.0f' % y for y in yticks])
            cb = plt.colorbar(im, ax=ax)
        
            if None:
                """ Quiver ion drift plot """
                x, y = pol2cart_vec(rad, theta, vi_N[::-1, :], -vi_E[::-1, :])
                ax.quiver(theta, rad, x, y)
                
                """ Quiver DMSP ion drift plot """
                rad = np.deg2rad(90) - np.deg2rad(dmsp[sat]['mlat'][tidx])
                theta = np.deg2rad(dmsp[sat]['mlong'][tidx])
                vi_N, vi_E = brng_mag_to_N_E(
                    np.deg2rad(dmsp[sat]['vi_dirn_MAG'][tidx]), dmsp[sat]['vi_mag'][tidx],
                )
                x, y = pol2cart_vec(rad, theta, -vi_N, vi_E)
                ax.quiver(theta, rad, x, y, color='m')

            """ Plot local noon dot """
            noon_glon = local_noon(time)  
            noon_mlon = glon_to_mlon(
                pot['Geographic Longitude'][-1, :], 
                pot['MLON (AACGM)'][-1, :], noon_glon,
            )
            noon_rad, noon_theta = get_rad_theta(latlim + 0.5, noon_mlon)

            ax.plot(noon_theta, noon_rad, '.r', markersize=20)

            ax.set_rmax(np.deg2rad(90 - latlim))
            plt.savefig(plot_fname)
            print('******************************* Saved to %s' % plot_fname)


        time += timestep 


""" Load funcs """
            
def load_dmsp(sat, stime, etime, dl_dmsp=False, ):
    """ Load DMSP data using pysat. Optionally download """


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

    return dmsp


def load_ampere(fname, dati=None):
    """ 
    load AMPERE data
     if no date/time requested assume that we're doing the whole
     AMPERE netcdf file
     otherwise, do just the requested time
    """
    dset = ncread_vars(fname)
    blocks = np.arange(dset['Jr'].shape[0])

    nlon = dset['nlon'][0] + 0 
    nlat = dset['nlat'][0] + 0 
    vars = {}
    for vn in ['dBnorth2', 'dBeast1', 'Jr']:
        vars[vn] = []
        for block in blocks:
            tmp = np.reshape(dset[vn][block], (nlon, nlat))
    
            # add pole:
            tmp = np.hstack((np.ones((tmp.shape[0], 1)) * tmp[:, 0].mean(), tmp)) 
            vars[vn].append(tmp.T)

    return vars


def load_igrf(lat, lon, time, refh=800):
    # Return IGRF coefficients in Tesla
    Be, Bn, Bu = igrf(lon, lat, refh, time)
    Be, Bn, Bu = Be * 1e-9, Bn * 1e-9, Bu * 1e-9 # nT -> T
    B0 = np.squeeze(np.sqrt(Be**2 + Bn**2 + Bu**2))
    Bu = np.squeeze(Bu)

    return B0, Be, Bn, Bu


""" Calculations/Processing """

def proc_dmsp(dmsp, dec_rate=1, mlat_cutoff=60, np_latlon=None):
    """ Process DMSP """
    assert mlat_cutoff > 0, 'SH filtering not implemented yet'
    try:
        dmsp_vals = dmsp.data[::dec_rate] # decimate
    except:
        return []
    dmsp_vals = dmsp_vals[(dmsp_vals['mlat'] > mlat_cutoff) & np.isfinite(dmsp_vals['hor_ion_v'])]
    dmsp_times = dmsp_vals.index
    lats, lons, alts = dmsp_vals['gdlat'], dmsp_vals['glon'], dmsp_vals['gdalt']
    brng = calc_bearings(lats[:-1], lons[:-1], alts[:-1], lats[1:], lons[1:])
    dmsp_vals = dmsp_vals[:-1]
    dmsp_vals['vi_mag'] = dmsp_vals['hor_ion_v']
    dmsp_vals['vi_dirn_geo'] = calc_vel_dirn(lats, lons, dmsp_times, brng, dmsp_vals['vi_mag'])
    dmsp_vals['vi_mag_model'] = np.ones(len(dmsp_vals)) * np.nan

    # add conversion to MAG drift directions if north pole is provided
    if np_latlon:
        np_bearings = calc_bearings(lats, lons, alts, 
            np.ones(lats.shape) * np_latlon[0], np.ones(lats.shape) * np_latlon[1])
        dmsp_vals['vi_dirn_MAG'] = dmsp_vals['vi_dirn_geo'] + np_bearings[:-1]

    return dmsp_vals


def calc_pflux(
    db_fname='data/ampere/20140522.0000.86400.600.north.grd.ncdf',
    pot_fname='data/pot_sami_cond/may14_euvac/ampere_mix_2014-05-22T00-10-00Z.nc',
    tind = 0
):
    """ !!!In progress!!! calculate Poynting flux from MIX """
    dB = load_ampere(db_fname)
    En, Ee, glat, glon, mlat, mlon = calc_mix_efld(pot_fname)


    dBn, dBe = [dB['dBnorth2'][tind], dB['dBeast1'][tind]]
    Pflux = En * dBe + Ee * dBn

    plt_vals = {'FAC': dB['Jr'][tind], 'Pot': pot['Potential'][:, :-1], 'Pflux': Pflux}
    lat, lon = pot['MLAT (AACGM)'][:, 0].copy(), pot['MLON (AACGM)'][0, :].copy()
    lon = lon[:-1].copy()

    fig, ax = plt.subplots(1, 1,subplot_kw={'projection': 'polar'})
    rad = np.deg2rad(90) - np.abs(np.deg2rad(lat))
    lon[lon < 0] = lon[lon < 0] + 360
    theta = np.deg2rad(lon)
    pv = dB['Jr'][0]

    # contour data over the map.
    im = ax.contourf(theta, rad, pv, vmin=-1, vmax=1)
    cb = plt.colorbar(im, ax=ax)
    plt.show()


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


def calc_bearings(latas, lonas, alts, latbs, lonbs):
    """ Calculate bearings using reference ellipsoid (deg/km) """
    wgs84 = nv.FrameE(name='WGS84')

    assert len(latas) == len(latbs), 'Assuming latas and latbs are equal-length arrays'

    brng_deg = np.zeros(len(latas)) * np.nan
    for ind, latb in enumerate(latbs):
        depth = -alts[ind] * 1E3
        pointB = wgs84.GeoPoint(latitude=latas[ind], longitude=lonas[ind], z=depth, degrees=True)
        lon = lonbs[ind]
        pointA = wgs84.GeoPoint(latitude=latb, longitude=lon, z=depth, degrees=True)
        p_AB_N = pointA.delta_to(pointB)  # note we want the bearing at point A
        brng_deg[ind] = p_AB_N.azimuth_deg 

    return brng_deg


def calc_mix_efld(pot):
    """ calculate the E-field from the MIX potential (take gradient in geo lat/lon)
    latE, lonE are geographic
    We also need AACGM lat/lon for intercomparison with other datasets
    """
    mix = {}
    pot_w = wrap(pot)

    mix['En'], mix['Ee'], mix['glat'], mix['glon'] = calc_efield(
        pot_w['Potential'] * 1E3,  # kV -> V
        pot_w['Geographic Latitude'],
        pot_w['Geographic Longitude'],
    )
    mix['mlat'] = pot['MLAT (AACGM)'][:, :-1]
    mix['mlon'] = pot['MLON (AACGM)'][:, :-1]

    return mix


def calc_exb_drift(efld, B0, Bu, lat, lon, s=0.1):
    """     
    Calculate the horizontal 'ExB' drift component in the specified direction 
    at the specified LAT/LON (geographic). 
    efld defined as in calc_mix_efld above
    Note we're assuming B is vertical
    """
    # TODO: Consider using the spline to calculate the derivatives directly

    # convert to radians, flatten and discard low-latitude boundary NaNs
    theta = np.pi - np.deg2rad(efld['glat'][1:,:]).flatten()
    phi = np.deg2rad(efld['glon'][1:,:]).flatten()
    data1 = efld['En'][1:,:].flatten()
    data2 = efld['Ee'][1:,:].flatten()

    theta_i = np.pi - np.deg2rad(lat)
    phi_i = np.deg2rad(lon)
    phi[phi < 0] += 2 * np.pi
    phi_i[phi_i < 0] += 2 * np.pi

    # interpolate
    intobj1 = scipy.interpolate.SmoothSphereBivariateSpline(theta, phi, data1, s=s)
    intobj2 = scipy.interpolate.SmoothSphereBivariateSpline(theta, phi, data2, s=s)

    En = intobj1(theta_i, phi_i, grid=False)
    Ee = intobj2(theta_i, phi_i, grid=False)

    # calculate ion drift vector 
    vi_E = En * Bu / B0**2
    vi_N = -Ee * Bu / B0**2

    vi = np.squeeze(np.sqrt(vi_E ** 2 + vi_N ** 2))
    vi_dir = np.squeeze(np.rad2deg(np.arctan(vi_E / vi_N)))

    return vi, vi_dir, vi_E, vi_N


def calc_efield(pot, glat, glon):
    """ 
    Calculate E-field from potential pot, using partial derivative approach 
    Assume pot is in Volts, output will be Volts/metre

    glat, glon used to calculate great-circle distances
    En, Ee defined in AACGM north/east direction
    """

    shapeE = (pot.shape[0] - 2, pot.shape[1] - 2)
    En = np.zeros(shapeE)
    Ee = np.zeros(shapeE)
    latE = np.zeros(shapeE)
    lonE = np.zeros(shapeE)
    wgs84 = nv.FrameE(name='WGS84')
    z = -780E3  # depth in m

    for lati in range(shapeE[0]):
        for loni in range(shapeE[1]):
            latE[lati, loni] = glat[lati+1, loni+1]
            lonE[lati, loni] = glon[lati+1, loni+1]

            pointA = wgs84.GeoPoint(
                latitude=glat[lati, loni+1],
                longitude=glon[lati, loni+1],
                z=z, degrees=True,
            )
            pointB = wgs84.GeoPoint(
                latitude=glat[lati+2, loni+1],
                longitude=glon[lati+2, loni+1],
                z=z, degrees=True,
            )
            pointC = wgs84.GeoPoint(
                latitude=glat[lati+1, loni],
                longitude=glon[lati+1, loni],
                z=z, degrees=True,
            )
            pointD = wgs84.GeoPoint(
                latitude=glat[lati+1, loni+2],
                longitude=glon[lati+1, loni+2],
                z=z, degrees=True,
            )

            gc_distLat_m = pointA.delta_to(pointB).length
            gc_distLon_m = pointC.delta_to(pointD).length

            En[lati, loni] = (pot[lati, loni + 1] - pot[lati + 2, loni + 1]) / gc_distLat_m  # direction of increasing AACGM lat
            Ee[lati, loni] = (pot[lati + 1, loni + 2] - pot[lati + 1, loni]) / gc_distLon_m  # direction of eastward AAGCM lon

    return En, Ee, latE, lonE


""" Helpers """

def brng_mag_to_N_E(vi_dir, vi_mag):
    """ go from vector magnitude and direction to N, E components """
    vi_N = np.cos(vi_dir) * vi_mag
    vi_E = np.sin(vi_dir) * vi_mag

    return vi_N, vi_E


def zero_360(vals):
    """ set values between 0 - 360 degrees """

    if isinstance(vals, float):
        vals = np.array(vals)
    vals[vals >= 360] -= 360
    vals[vals < 0] += 360
    return vals


def pol2cart_vec(rad, theta, dr, dt):
    """ matplotlib polar quiver not implemented properly... """
    x = dr * np.cos(theta) - dt * np.sin (theta)
    y = dr * np.sin(theta) + dt * np.cos(theta)

    return(x, y)


def wrap(dictarr):
    """ 
    Wrap grid in longitude and above pole/below lowest lat for gradient calculation
    """
    d2 = {}
    for k, arr in dictarr.items():
        # add opposite longitude ring above the pole
        half_lon = int(arr.shape[1] / 2)
        lon_idx = np.concatenate([np.arange(half_lon, arr.shape[1]), np.arange(half_lon)])
        arr2 = np.concatenate([arr[1, [lon_idx,]], arr], axis=0)

        # add duplicate ring below the lowest latitude
        arr3 = np.concatenate([arr2, arr2[[-1,], :]], axis=0)

        # Add another longitude wrap before the first one
        arr4 = np.concatenate([arr3[:, [-2,]], arr3], axis=1)
        d2[k] = arr4

    return d2


def grad_sph(sfld, theta, phi, rad=65E3):
    """ !!!in progress!!! take gradient of scalar field in theta/phi directions """ 
    d_sfld = np.gradient(sfld, [theta, phi])
    d_sfld_dtheta = 1 / r * d_sfld[0]
    d_sfld_dphi = 1 / (r * np.sin(theta)) * d_sfld[1]
    
    
def local_noon(time):
    """ Calculate local noon longitude from UTC time """
    noon_lon = zero_360(-(time.hour / 24 + time.minute / 1440) * 360 + 180)

    return noon_lon
    

def glon_to_mlon(glon, mlon, glon_i):
    """ interpolate geographic longitude to AACGM lon """
    intf = scipy.interpolate.interp1d(zero_360(glon), zero_360(mlon), fill_value='extrapolate') 
    noon_mlon = intf(zero_360(glon_i))

    return noon_mlon


def get_rad_theta(lat, lon):
    """ convert lat, lon (deg) to r/theta for polar plotting """
    rad = np.deg2rad(90) - np.abs(np.deg2rad(lat))
    theta = np.deg2rad(lon)

    return rad, theta 


if __name__ == '__main__':
    main()




















