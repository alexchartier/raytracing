import numpy as np
import time
import os
import matplotlib.pyplot as plt
import h5py
import datetime as dt
import nvector as nv
import pandas as pd
import pymap3d as pm
import nc_utils

wgs84_pm = pm.Ellipsoid('wgs84')
wgs84 = nv.FrameE(name='WGS84')
RAD_EARTH = 6371E3

""" Utilities to deal with MIT Haystack GPS/GNSS files. Read functions come from Bill Rideout PDF (MIT Haystack) """

def load_gps(gps_fn, stime, timestep=dt.timedelta(minutes=60)):
    """ Read the MIT Haystack GPS file and downsample """
    gps_fn = os.path.abspath(os.path.expanduser(gps_fn))
    assert stime.tzinfo == dt.timezone.utc, 'Provide UTC times'

    # Define the processing times
    etime = stime + dt.timedelta(days=1)
    times = []
    t = stime
    while t < etime:
        times.append(t)
        t += timestep

    # Get the times: We only need to do this step once, so that getting data from other sites is faster
    print('Initial GPS load (takes 2 - 3 mins): %s' % gps_fn)
    with h5py.File(gps_fn, 'r') as f:
        timeArr = f['Data']['Table Layout']['ut1_unix']
    print('... GPS loaded')

    # get the list of unique times from this array
    gps_timestamps = np.unique(timeArr)
    timestamps = [t.timestamp() for t in times]
    
    # Run through the requested times and store out data
    gps_in = {}
    for t in timestamps:
        print(t)
        assert t in gps_timestamps, 'requested time not in file - check'
        data_t = getData_t(gps_fn, timeArr, t)
        vn = [v for v in data_t.dtype.fields.keys()]

        if t == timestamps[0]:
            for v in vn:
                gps_in[v] = []
        for entry in data_t:
            for ind, v in enumerate(vn):            
                gps_in[v].append(entry[ind])

    gps = {}
    gps['time'] = []
    for ti, year in enumerate(gps_in['year']):
        gps['time'].append(dt.datetime(
            year, gps_in['month'][ti], gps_in['day'][ti], gps_in['hour'][ti],
            gps_in['min'][ti], gps_in['sec'][ti],
        ))

    gps['time'] = np.array(gps['time'])

    wanted_vn = 'gps_site', 'sat_id', 'gdlatr', 'gdlonr', 'los_tec', 'dlos_tec', 'azm', 'elm'
    for v in wanted_vn:
        gps[v] = np.array(gps_in[v])

    gps = pd.DataFrame.from_dict(gps)
    gps.set_index('time')

    return gps


def calc_tx_rx_coords(gps, ref_rx_alt_m=0.):
    """ Determine the transmitter and receiver XYZ coords based on receivers at 0 m alt """
    rx_XYZ_full = []
    tx_XYZ_full = []
    for idx in range(len(gps)):
        print('Processing %i of %i' % (idx + 1, len(gps)))
        gps_i = gps.iloc[idx]
        rx_XYZ = np.squeeze(sphcart(ref_rx_alt_m, gps_i['gdlatr'], gps_i['gdlonr']))
        tx_XYZ = get_tx_loc(rx_XYZ, gps_i['azm'], gps_i['elm'], gps_i['gdlatr'], gps_i['gdlonr']) 
        rx_XYZ_full.append(rx_XYZ)
        tx_XYZ_full.append(tx_XYZ)
    
    gps['rx_XYZ'] = rx_XYZ_full
    gps['tx_XYZ'] = tx_XYZ_full

    return gps


def get_tx_loc(rx_XYZ, az, el, lat, lon, ref_rx_alt_m=0., ref_tx_alt_m=2E7):
    """ Calculate transmitter location from receiver location and az/el"""

    # get the 1000-km unit vector towards the satellite
    unitvec_XYZ_loc = np.array(pm.aer2ecef(az, el, 1E6, lat, lon, ref_rx_alt_m, wgs84_pm)) - rx_XYZ

    # Calculate the intersection of that vector with a 20k-km altitude sphere
    tx_XYZ = sphere_intersect(RAD_EARTH + ref_tx_alt_m, unitvec_XYZ_loc, rx_XYZ)

    return tx_XYZ


def sphere_intersect(rad, line, origin, c=0):
    """
    Location of intersection of a line with a sphere
     rad - radius of sphere
     line - vector of line
     origin - origin of line
     c - centre of sphere
     dist - distance along unit vector to point

    Should return the closest location that intersects
    """

    mod_l_2 = np.multiply(line, line).sum(0)  # Dot product 
    oc = origin - c 
    mod_oc_2 = np.multiply(oc, oc).sum(0)
    t1 = -2 * np.multiply(line, oc).sum(0)
    t2 = np.sqrt((2 * np.multiply(line, oc).sum(0)) ** 2 - 4 * mod_l_2 * (mod_oc_2 - rad ** 2)) 
    d1 = (t1 + t2) / (2 * mod_l_2)
    d2 = (t1 - t2) / (2 * mod_l_2)

    assert (d1 > 0).all(), 'Some d1s negative - calculation must have failed'
    pt = origin + line * d1

    return pt


def getData_t(madrigalFile, timeArr, unixTime):
    """getData_t returns a numpy recarray of all the data in
    madrigalFile from a given time

    Inputs:
    madrigalFile - path to Madrigal LOS Hdf5 file
    timeArr - a numpy array of the times (ut1_unix) in the entire
    file.
    unixTime - the unix time to select

    Returns a numpy recarray of the subset of data from the file from
    that one site
    """
    with h5py.File(madrigalFile, 'r') as f:
        # get a list of all the indices with the right site
        time_indices = timeArr == unixTime
        timeData = f['Data']['Table Layout'][time_indices]

    return(timeData)


def getGpsData_t(madrigalFile, timeArr, satTypeArr, unixTime):
    """getTimeGpsData returns a numpy recarray of all the data in
    madrigalFile from a given time
    from the GPS satellite constellation only

    Inputs:
    madrigalFile - path to Madrigal LOS Hdf5 file
    timeArr - a numpy array of the times (ut1_unix) in the entire
    file.
    atTypeArr - a numpy.array of 8 char strings 'GPS ' or
    'GLONASS ' (others possible in future)
    unixTime - the unix time to select

    Returns a numpy recarray of the subset of data from the file from
    that one site
    """
    with h5py.File(os.path.expanduser(madrigalFile), 'r') as f:
        # get a list of all the indices with the right site
        # use b'GLONASS ' to get GLONASS only
        indices = np.logical_and(timeArr == unixTime, satTypeArr == b'GPS ')
        timeData = f['Data']['Table Layout'][indices]
    return(timeData)


def test_time_read(gps_fn):
    """ Test out the time read speed to demonstrate initial read is the only slow one """
    
    gps_fn = os.path.expanduser(gps_fn)
    t = time.time()
    # we only need to do this step once, so that getting data from other sites is faster
    f = h5py.File(gps_fn, 'r')
    timeArr = f['Data']['Table Layout']['ut1_unix']

    # get the list of unique times from this array
    times = np.unique(timeArr)
    f.close()
    timeData = getData_t(gps_fn, timeArr, times[0])
    print('Took %f secs to get %i measurements from the first time' % (time.time()-t, len(timeData)))

    # show that accessing the second time is faster
    t = time.time()
    timeData = getData_t(gps_fn, timeArr, times[1])
    print('Took %f secs to get %i measurements from the second time' % (time.time()-t, len(timeData)))
    print('The following are the column names and data types in this file:')
    for colName, colType in timeData.dtype.fields.items():
        print('%s\t%s' % (colName, colType))
    print('For example, here are the first four rows:')
    for i in range(4):
        print(i, timeData[i])
    print('For example, here are the first four rows of line of site TEC:')
    for i in range(4):
        print(i, timeData['los_tec'][i])
    print('There are %i unique times in the file' % (len(times)))


def test_gps_read(gps_fn):
    """ Test out the GPS read speed to demonstrate initial read is the only slow one """
    gps_fn = os.path.expanduser(gps_fn)
    t = time.time()

    # we only need to do this step once, so that getting data from other sites is faster
    f = h5py.File(gps_fn, 'r')
    timeArr = f['Data']['Table Layout']['ut1_unix']

    # get the list of unique times from this array
    times = np.unique(timeArr)
    try:
        satTypeArr = f['Data']['Table Layout']['gnss_type']
    except:
        # if this column does not exist, ALL data is GPS (gnss_type == 0)
        satTypeArr = np.zeros((len(timeArr),), dtype='S8')
        satTypeArr[:] = 'GPS '

    f.close()
    timeData = getGpsData_t(gps_fn, timeArr, satTypeArr, times[0])
    print('Took %f secs to get %i measurements from the first time' % (time.time()-t, len(timeData)))

    # show that accessing the second time is faster
    t = time.time()
    timeData = getGpsData_t(gps_fn, timeArr, satTypeArr, times[1])
    print('Took %f secs to get %i measurements from the second time' % (time.time()-t, len(timeData)))
    print('There are %i unique times in the file' % (len(times)))


def gen_sitelist(gps, npts_wanted=100):
    """ Grab all the unique sites in GPS df and pull out npts of them, based on max separation """
    slist_full = get_sitelist(gps)
    all_lats = slist_full['gdlatr'].tolist()
    all_lons = slist_full['gdlonr'].tolist()
    all_sites = slist_full['sites'].tolist()
    good_lats = [all_lats.pop(0),]
    good_lons = [all_lons.pop(0),]
    good_sites = [all_sites.pop(0),]
    
    while len(good_lats) < npts_wanted:
        distmat = np.ones((len(good_lats), len(all_lats)))
        for i1, lat1 in enumerate(good_lats):
            distmat[i1, :] = calc_greatcircle_dists(lat1, good_lons[i1], all_lats, all_lons) 

        dist = np.prod(distmat / distmat.max(), axis=0)
        idx = np.argmax(dist)
        good_lats.append(all_lats.pop(idx))
        good_lons.append(all_lons.pop(idx))
        good_sites.append(all_sites.pop(idx))

    slist = {}
    slist['sites'] = np.array(good_sites)
    slist['gdlatr'] = np.array(good_lats)
    slist['gdlonr'] = np.array(good_lons)
    slist['XYZ'] = sphcart(0, np.array(good_lats), np.array(good_lons))
    slist['LLA'] = sphcart(0, np.array(good_lats), np.array(good_lons))

    #plt.plot(slist['gdlonr'], slist['gdlatr'], '.k'); plt.plot(good_lons, good_lats, '.r'); plt.show()
    return slist


def calc_greatcircle_dists(lat1, lon1, lat2, lon2):
    """ Great circle distances based on WGS84 ellipsoid """
    pointA = wgs84.GeoPoint(latitude=lat1, longitude=lon1, z=0, degrees=True)
    pointB = wgs84.GeoPoint(latitude=lat2, longitude=lon2, z=0, degrees=True)
    dist = pointA.delta_to(pointB).length
    
    return dist


def sphcart(alt, lat, lon, degrees=True):
    """ spherical (km, deg) to cartesian (m) """
    cart = wgs84.GeoPoint(
        latitude=lat, longitude=lon, z=-alt * 1E3, degrees=degrees,
    ).to_ecef_vector().pvector

    return cart


def cartsph(cart, degrees=True):
    """ cartesian (m) to spherical (km, deg) """
    n_EB_E, z_EB = nv.p_EB_E2n_EB_E(cart)
    alt = -z_EB / 1E3 
    lat_EB, lon_EB = nv.n_E2lat_lon(n_EB_E)
    if degrees:
        lat = np.rad2deg(lat_EB)
        lon = np.rad2deg(lon_EB)

    return alt, lat, lon


def get_sitelist(gps):
    """ List of unique sites from GPS df """
    slist = {}
    [slist['sites'], idx] = np.unique(gps['gps_site'], return_index=True)
    for k in 'gdlatr', 'gdlonr':
        slist[k] = gps[k][idx] 

    return slist


def test_conversions(tol=1E-3):
    """ check cartsph and sphcart """
    alt = 6400.
    lat = 20. 
    lon = -20.

    cart = sphcart(alt, lat, lon)
    [alt_out, lat_out, lon_out] = cartsph(cart)

    assert np.abs(alt - alt_out) < tol, 'alt wrong'
    assert np.abs(lat - lat_out) < tol, 'lat wrong'
    assert np.abs(lon - lon_out) < tol, 'lon wrong'


def downsample_to_slist(gps, slist):
    """ Throw out the sites not in slist """
    goodidx = gps['gps_site'].values == slist['sites'][0]
    for site in slist['sites']:
        goodidx += gps['gps_site'].values == site

    gps_short = gps[goodidx]
    print('Downsampled the GPS to include only the sitelist entries - %i' % len(gps_short))
    return gps_short


def cleanup(gps):
    """ TODO e.g. remove negatives, maybe index by time, specify min nsats=4 """
    neg_vals = gps['los_tec'] < 0
    gps = gps.drop(gps.index[neg_vals])
    return gps


def load_and_preproc_mit_gps(gps_fn, slist_pkl_fn, time):
    """ Run through from MIT GPS input file to short, clean output """
    gps = load_gps(gps_fn, time)
    slist = nc_utils.unpickle(slist_pkl_fn)
    gps_short = downsample_to_slist(gps, slist)
    gps_short_XYZ = calc_tx_rx_coords(gps_short)
    return cleanup(gps_short_XYZ)


def write_sitelists(gps, slist, slist_fn, slist_fn_short):
    
    # write the full sitelist
    idx = np.unique(gps['gps_site'], return_index='True')
    unique_gps = gps.iloc[idx[1]]
    unique_gps = unique_gps.drop(['time', 'sat_id', 'los_tec', 'dlos_tec', 'azm', 'elm'], axis=1)
    unique_gps_str = unique_gps.to_string(header=False, index=False)
    with open(slist_fn, 'w') as f:
        f.write(unique_gps_str)


    # write the short sitelist
    del slist['XYZ']
    del slist['LLA']
    gps_short = pd.DataFrame(slist)
    gps_str = gps_short.to_string(header=False, index=False)
    with open(slist_fn_short, 'w') as f:
        f.write(gps_str)


if __name__ == '__main__':
    gps_fn_fmt = '~/data/gps/mit_hdf/los_%Y%m%d.001.h5'
    slist_pkl_fn = '~/data/gps/sitelists/global_150.pkl'

    stime = dt.datetime(2019, 3, 1, tzinfo=dt.timezone.utc)
    etime = dt.datetime(2019, 3, 2, tzinfo=dt.timezone.utc)

    """
    # Get the GPS preprocessed
    time  = stime
    while time < etime:
        gps_fn = time.strftime(gps_fn_fmt)
        gps_pkl_fn = time.strftime(gps_pkl_fn_fmt)
        gps = load_gps(gps_fn, time)
        nc_utils.pickle(gps, gps_pkl_fn)
        time += dt.timedelta(days=1)

    """ #Load pkls instead

    gps_pkl_fn_fmt = '~/data/gps/pkl/los_%Y%m%d.001.pkl'
    gps_pkl_fn_fmt_2 = '~/data/gps/pkl/los_small_%Y%m%d.001.pkl'
    gps_pkl_fn_fmt_3 = '~/data/gps/pkl/los_small_XYZ_%Y%m%d.001.pkl'
    gps_pkl_fn = stime.strftime(gps_pkl_fn_fmt)
    gps_pkl_fn_2 = stime.strftime(gps_pkl_fn_fmt_2)
    gps_pkl_fn_3 = stime.strftime(gps_pkl_fn_fmt_3)
    gps = nc_utils.unpickle(gps_pkl_fn)
    gps_short = nc_utils.unpickle(gps_pkl_fn_2)
    slist = nc_utils.unpickle(slist_pkl_fn)

    # write some info for matlab to plot
    write_sitelists(gps, slist, 'data/gps/sitelist.txt', 'data/gps/sitelist_150.txt')

    # Calculate the sitelist
    slist = gen_sitelist(gps, 150)
    nc_utils.pickle(slist, slist_pkl_fn)
    

    # Downsample the GPS to the sitelist, and calculate coordinates 
    gps_short = downsample_to_slist(gps, slist)
    nc_utils.pickle(gps, gps_pkl_fn_2)
   
    # Calculate the XYZ coordinates 
    gps_short_XYZ = calc_tx_rx_coords(gps_short)

    # clean up the data
    gps_short_XYZ = cleanup(gps_short_XYZ)
    
    nc_utils.pickle(gps_short_XYZ, gps_pkl_fn_3)












