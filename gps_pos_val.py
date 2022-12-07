"""
gps_pos_val.py
Calculate the GPS ionospheric positioning error based on files from MIT Haystack
"""
import numpy as np
import datetime as dt
import os
import pickle
from nc_utils import ncread_vars, load_nc
import scipy.spatial.qhull as qhull
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import nvector as nv
import matplotlib.dates as mdates
import h5py
wgs84 = nv.FrameE(name='WGS84')
xfmt = mdates.DateFormatter('%m/%d %H:%M')


def main(sami_fn, gps_fn, out_fn):

    sami = load_sami(sami_fn)

    #ds = downsample_hdf5_file(gps_fn)

    try: 
        gps = load_preprocessed_gps(preproc_gps_fn)
    except:
        gps_full = load_gps(gps_fn)
        
        breakpoint()
        sites = gen_sitelist(gps_full)  # just set it to the 1st one first
        
        preproc_gps(gps, sitelist)


    # Run through selected sites and calculate position errors
    rx = {}
    for site_name, site_fname in sites.items():
        xyz = np.array([sitef.receiver_x * 1E3, sitef.receiver_y * 1E3, sitef.receiver_z * 1E3])
        alt, lat, lon = cartsph(xyz)
        rx[site_name] = {'lla': [lat.tolist(), lon.tolist(), alt.tolist()], 'XYZ': xyz} 

        # data
        data = proc_obs_delays(times, site_fname)  
        data = calc_pos_err(data)

        # model
        modvals = {}
        in_fn = os.path.join(v, 'prior_reg/%Y%j-%H%M_prior_reg.nc')
        rxi = rx[site_name]
        rxi['Lat'], rxi['Lon'], rxi['Alt'] = rxi['lla']
        mod = proc_mod_delays(
            times, data['tx'], rxi, data['prn'],
            in_fname=in_fn, 
            site_name=site_name, 
        )


def proc_obs_delays(times, site_fname):
    # Calculate the observed ionospheric positioning error based on observed TEC
    out = {
        'delay': [],
        'sTEC': [],
        'tx': [],
        'rx': [],
        'prn': [],
        'elv': [],
        'times': [],
    }
    for tind, time in enumerate(times):
        # Load slant TEC each day
        if ((time.hour == 0) and (time.minute == 0)) or tind == 0:
            # Load the data at the beginning of the day
            day = time
            site_fname_t = time.strftime(site_fname)
            data = ncread_vars(site_fname_t)
            params = load_nc(site_fname_t)
            rx_XYZ = np.array([params.receiver_x, params.receiver_y, params.receiver_z]) * 1E3
            rx_LLA = cartsph(rx_XYZ)
            for k, v in data.items():
                data[k] = np.array(v)
            secs = (data['times'] - np.floor(data['times'][0])) * 86400
            if secs.max() > 86400:  # more than a whole day of data - assume it started on previous day
                secs -= 86400
            tx_XYZ = np.array([data['x_position'], data['y_position'], data['z_position']]).T * 1E3
            delays = data['TECData'] / 6.13
            if data['TECData'].min() < -25:
                pdb.set_trace()

        tind = []
        absd = np.abs(secs - (time - day).total_seconds())
        tind = absd == absd.min()
        if absd.min() < 60:  # 'Must be within a minute'
            out['sTEC'].append(data['TECData'][tind])
            out['delay'].append(delays[tind])
        else:
            out['sTEC'].append(data['TECData'][tind] * np.nan)
            out['delay'].append(delays[tind] * np.nan)

        out['prn'].append(data['satellites'][tind])
        tx_XYZ_t = tx_XYZ[tind, :]
        out['tx'].append({'XYZ': tx_XYZ_t})
        out['rx'].append({'XYZ': rx_XYZ})
        elv = np.zeros(sum(tind)) * np.nan
        for ind in range(len(elv)):
            elv[ind] = pymap3d.ecef2aer(*tx_XYZ_t[ind, :], *rx_LLA)[1]
        if len(elv) != len(delays[tind]):
            pdb.set_trace()
        out['elv'].append(elv)

    out['times'] = times

    return out


def proc_mod_delays(
    times, tx, rx, prn,
    in_fname='/Users/chartat1/fusionpp_data/sami/prior_reg/%Y%j-%H%M_prior_reg.nc',
    site_name='aoml',
):
    # Calculate slant TEC and delays from the model 
    out = {
        'delay': [],
        'sTEC': [],
        'tx': [],
        'rx': [],
        'prn': [],
        'elv': [],
        'azm': [],
        'times': [],
    }
    for tind, time in enumerate(times): 
        print('Pre-processing ...')
        print(time.strftime('%Y/%b/%d %H:%M'))

        # Load model and Delaunay grid
        mod = ncread_vars(time.strftime(in_fname))
        if tind == 0:
            alt, lat, lon = np.meshgrid(mod['alt'], mod['lat'], mod['lon'], indexing='ij')
            mod_XYZ = sphcart(alt.flatten(), lat.flatten(), lon.flatten())
            mod_tri = qhull.Delaunay(mod_XYZ.T)
        mod['XYZ'] = mod_XYZ
        mod['tri'] = mod_tri
        mod['dene'] *= 1E6   # cm3 to m3

        # Calculate the delay
        azm, elv = calc_az_el(prn[tind], tx[tind]['XYZ'], rx['lla'])
        sTEC = calc_sTEC(tx[tind]['XYZ'], rx['XYZ'], mod) / 1E16
        delays = sTEC / 6.13 

        # Store data
        out['delay'].append(delays)
        out['sTEC'].append(sTEC)
        out['tx'].append(tx[tind])
        out['rx'].append(rx)
        out['prn'].append(prn[tind])
        out['elv'].append(elv)
        out['azm'].append(azm)
        out['times'].append(time)

    out['times'] = np.array(out['times']) 
    return out


def calc_pos_err(out):
    # Calculate positioning errors 
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


def get_pos_err(tx, rx, elv, delays):
    # Calculate positioning error
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
    res = minimize(cost_func, rx, method='Powell', options={'ftol':0.01})
    pos_err = np.sqrt(np.sum((rx - res.x) ** 2))
  
    print('delays')
    print(*['%1.1f m, ' % d for d in delays])
    if pos_err > 50:
        pdb.set_trace()
    print('3D position error: %2.1f m' % pos_err)
    print('Inversion status: %s' % res.message)

    return pos_err


def calc_sTEC(tx, rx, mod):
    # Get the delay on each PRN
    # tx: nx3 array of ECEF transmitter coords 
    # rx: nx3 array of ECEF receiver coords 
    npts = tx.shape[0]
    sTEC = np.zeros((npts,))
    for ind in range(npts):
        sTEC[ind] = get_slant_TEC(mod, tx[ind, :].copy(), rx.copy())
    return sTEC


def sTEC_to_L1_delay(sTEC):
    delay = sTEC / (6.13 * 1E16)
    return delay


def get_slant_TEC(mod, tx, rx):
    # Get the TEC from the model input file
    stepsize = 10E3
    xi = create_raypath(tx, rx, stepsize) 
    assert np.sqrt(np.sum(mod['XYZ'].T ** 2, 1)).max() > np.sqrt(np.sum(xi ** 2, 1)).max(), 'Ray goes above model top'
    ray_ne = interpolate(mod['dene'], mod['tri'], xi)
    sTEC = np.sum(ray_ne) * stepsize
    
    return sTEC


def interpolate(vals, tri, uvw):
    d = 3
    simplex = tri.find_simplex(uvw)
    vtx = np.take(tri.simplices, simplex, axis=0)
    temp = np.take(tri.transform, simplex, axis=0)
    delta = uvw - temp[:, d]
    bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)
    wts = np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))

    return np.einsum('nj,nj->n', np.take(vals, vtx), wts)


def create_raypath(tx_XYZ, rx_XYZ, stepsize=1E3, modtop=7600E3):
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


def get_tx_rx_pos(datacoll_fname, sitename):
    # Get the satellite and receiver positions from the data collector
    tx = {}
    rx = {}
    PRN = None
    elv = np.nan
    datacoll = load_datacollector(datacoll_fname)
    datacoll['GPS']['sitename'] = np.array([s.lower() for s in datacoll['GPS']['sitename']])
    siteind = datacoll['GPS']['sitename'] == sitename.lower()
    if np.sum(siteind) == 0:
        print('Site %s not found in %s' % (sitename, datacoll_fname))
        return tx, rx, PRN, elv

    crds = 'Lat', 'Lon', 'Alt'
    for crd in crds:
        tx[crd] = datacoll['GPS']['end%s' % crd][siteind].astype(np.float)
        rx[crd] = datacoll['GPS']['start%s' % crd][siteind].astype(np.float)
    PRN = datacoll['GPS']['PRN'][siteind].astype(np.float)

    tx['XYZ'] = sphcart(tx['Alt'], tx['Lat'], tx['Lon']).T
    rx['XYZ'] = sphcart(rx['Alt'], rx['Lat'], rx['Lon']).T

    assert len(np.unique(rx['XYZ'])) == 3, 'Assuming the Rx points are unique'
    rx['XYZ'] = rx['XYZ'][0, :]
    for crd in crds:
        rx[crd] = rx[crd][0]

    npts = tx['XYZ'].shape[0]
    elv = np.zeros(npts) * np.nan
    for ind in range(npts):
        elv[ind] = pymap3d.ecef2aer(*tx['XYZ'][ind, :], rx['Lat'], rx['Lon'], rx['Alt'])[1]

    return tx, rx, PRN, elv


def calc_az_el(PRN, tx_XYZ, rx_lla):
    azm = np.zeros(PRN.shape) * np.nan
    elv = np.zeros(PRN.shape)
    for ind, prn in enumerate(PRN):
        azm[ind], elv[ind] = pymap3d.ecef2aer(*tx_XYZ[ind, :], *rx_lla)[:2]
    return azm, elv


def get_tx_from_sp3(sp3_fname, rx_lla, time):
    tx, PRN = read_sp3(sp3_fname, time)
    azm, elv = calc_az_el(PRN, tx['XYZ'], rx_lla)
    elvind = elv > 0
    tx['XYZ'] = tx['XYZ'][elvind, :]
    PRN = PRN[elvind]  
    return tx, PRN, elv[elvind]


def read_sp3(sp3_fname, time):
    with open(sp3_fname, 'r') as f:
        txt = f.readlines()
  
    PRN = [] 
    tx = {'XYZ': []}
    for ind, ln in enumerate(txt):
        if ln[0] == '*':  # time lines start with *
            ln_t = [int(l) for l in ln.split()[1:6]]
            if time == dt.datetime(*ln_t):  # look for matching time
                ct = 1
                while ct < 32:
                    line = txt[ind + ct].split()
                    if (line[0] == '*') or (line[0] == 'EOF'): # we have reached the next time or end of file
                        break
                    if line[0][:2] == 'PG':  # The new format
                        PRN.append(int(line[0][2:]))
                        tx['XYZ'].append([float(l) * 1E3 for l in line[1:4]])
                    else:  # the old format
                        PRN.append(int(line[1]))
                        tx['XYZ'].append([float(l) * 1E3 for l in line[2:5]])

                    ct += 1
            
    tx['XYZ'] = np.array(tx['XYZ'])
    PRN = np.array(PRN)
    assert len(PRN) < 33, 'Only 32 GPS satellites flying'
    return tx, PRN


def sphcart(alt, lat, lon):
    # km, deg to XYZ in m
    cart = wgs84.GeoPoint(
        latitude=lat, longitude=lon, z=-alt * 1E3, degrees=True,
    ).to_ecef_vector().pvector

    return cart


def cartsph(cart):
    # XYZ (m) to LLA (km/deg)

    n_EB_E, z_EB = nv.p_EB_E2n_EB_E(cart)
    alt = -z_EB / 1E3
    lat_EB, lon_EB = nv.n_E2lat_lon(n_EB_E)
    lat = np.rad2deg(lat_EB)
    lon = np.rad2deg(lon_EB)
    
    return alt, lat, lon


def test_conversions(tol=1E-3):
    alt = 6400.
    lat = 20.
    lon = -20.

    cart = sphcart(alt, lat, lon)
    [alt_out, lat_out, lon_out] = cartsph(cart)

    assert np.abs(alt - alt_out) < tol, 'alt wrong'
    assert np.abs(lat - lat_out) < tol, 'lat wrong'
    assert np.abs(lon - lon_out) < tol, 'lon wrong'
    

def load_sami(sami_fn):
    sami = ncread_vars(sami_fn)
    return sami


def load_gps(gps_fn):
    # Read the MIT Haystack GPS file
    gps_fn = os.path.abspath(os.path.expanduser(gps_fn))
    hf = h5py.File(gps_fn, 'r')
    hlay = hf["Data/Table Layout"]
    keys = hf.keys()
    if 'gnss_type' in keys:
        gtype = hlay['gnss_type']
    else:
        gtype = None
    gps = {}  
    breakpoint()
    gps['hour'] = hlay['hour'] * 1.0 
    gps['min'] = hlay['min'] * 1.0 
    gps['sec'] = hlay['sec'] * 1.0 
    
    gps['oprns'] = hlay['sat_id']
    gps['orlat'] = hlay['gdlatr']
    gps['orlon'] = hlay['gdlonr']
    gps['otec']  = hlay['los_tec']
    gps['odtec'] = hlay['dlos_tec']
    gps['oazm']  = hlay['azm']
    gps['oelm']  = hlay['elm']
    gps['rsite']  = hlay['gps_site']
    gps['ursite'] = np.unique(gps['rsite'])
    gps['nsites'] = gps['ursite'].size
    
    print("Found %i unique sites" % gps['nsites'])
    return gps


def calc_sitelist(gps):
    return sitelist


def downsample_hdf5_file(in_fname, pts_wanted=1000):
    #Get length of files and prepare samples
    source_file = h5py.File(os.path.expanduser(in_fname), "r")
    dataset = source_file['/Data/Table Layout']
    indices = np.sort(np.random.choice(dataset.shape[0], int(pts_wanted), replace=False))

    #checking we're extracting a subsample
    if pts_wanted > dataset.shape[0]:
        raise ValueError("Can't extract more rows than dataset contains. Dataset has %s rows" % dataset.shape[0] )

    #target_file =  h5py.File(out_fname, "w")
    for k in source_file.keys():
        dataset = source_file[k]
        breakpoint()
        dataset = dataset[indices,:,:,:]
    #    dest_dataset = target_file.create_dataset(k, shape=(dataset.shape), dtype=np.float32)
    #dest_dataset.write_direct(dataset)
    #target_file.close()
    source_file.close()


if __name__ == '__main__':

    test_conversions()

    sami_fn = '~/data/sami3/2019/sami3_regulargrid_elec_density_2019Mar01.nc'
    gps_fn = '~/data/gps/los_20190301.001.h5'
    out_fn = 'out.nc'
    main(sami_fn, gps_fn, out_fn)



























