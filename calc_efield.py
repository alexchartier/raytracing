import matplotlib.pyplot as plt
import numpy as np
import nvector as nv
from nc_utils import ncread_vars
from scipy.io import netcdf
from ppigrf import igrf
import scipy.interpolate


def main(
    db_fname='data/ampere/20140522.0000.86400.600.north.grd.ncdf',
    pot_fname='data/pot_sami_cond/may14_euvac/ampere_mix_2014-05-22T00-10-00Z.nc',
):
    tind = 0
    dB = get_data(db_fname)
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

    """ 
    fig, ax = plt.subplots(len(plt_vals), 1,subplot_kw={'projection': 'polar'})
    ct = 0
    for k, v in plt_vals.items(): 
        print(k)
        plot_polar(ax[ct], np.deg2rad(lat), np.deg2rad(lon), v)
        ct += 1 
    plt.show()
    """

def get_exb_drift_2(pot_fname, B, lat, lon):
    """ Calculate the E-field at the interpolation locations. 
    Use spherical spline interpolation with second derivatives
    TODO: finish this function
    """
    pot = ncread_vars(pot_fname)

    theta = np.pi - np.deg2rad(pot['Geographic Latitude']).flatten()
    theta_i = np.pi - np.deg2rad(lat).flatten()
    phi = np.deg2rad(pot['Geographic Longitude']).flatten()
    phi_i = np.deg2rad(lon).flatten()
    phi[phi < 0] += 2 * np.pi
    data = pot['Potential'].flatten()
    intobj = scipy.interpolate.SmoothSphereBivariateSpline(theta, phi, data, s=0.5)
    En, Ee = intobj(theta_i, phi_i, dtheta=1, dphi=1)


def calc_mix_efld(pot_fname):
    """ calculate the E-field from the MIX potential (take gradient in geo lat/lon)
    latE, lonE are geographic
    We also need AACGM lat/lon for intercomparison with other datasets
    """
    mix = {}
    pot = ncread_vars(pot_fname)
    pot_w = wrap(pot)

    mix['En'], mix['Ee'], mix['glat'], mix['glon'] = calc_efield(
        pot_w['Potential'] * 1E3,  # kV -> V
        pot_w['Geographic Latitude'], 
        pot_w['Geographic Longitude'],
    )
    mix['mlat'] = pot['MLAT (AACGM)'][:, :-1]
    mix['mlon'] = pot['MLON (AACGM)'][:, :-1]

    return mix


def get_exb_drift(efld, B0, Bu, dirn, lat, lon, s=0.1):
    """     
    Calculate the horizontal 'ExB' drift component in the specified direction 
    at the specified LAT/LON (geographic). 
    efld defined as in calc_mix_efld above
    Note we have to convert AACGM to LAT/LON at the right altitude before using this!
    Also note we're assuming B is vertical
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
    breakpoint()
    En = intobj1(theta_i, phi_i, grid=False)        
    Ee = intobj2(theta_i, phi_i, grid=False)

    # calculate ion drift vector 
    vi_E = En * Bu / B0**2
    vi_N = -Ee * Bu / B0**2

    vi = np.sqrt(np.sum(vi_E ** 2 + vi_N ** 2))
    vi_dir = np.rad2deg(np.arctan(vi_E / vi_N))

    return vi, vi_dir
  

def get_igrf(lat, lon, time, refh=800):
    # Return IGRF coefficients in Tesla
    Be, Bn, Bu = igrf(lon, lat, refh, time)
    Be, Bn, Bu = Be * 1e-9, Bn * 1e-9, Bu * 1e-9 # nT -> T
    B0 = np.sqrt(Be**2 + Bn**2 + Bu**2).reshape((1, -1))
    Bu = Bu.reshape((1, -1))

    return B0, Be, Bn, Bu


def plot_polar(
        ax, lat, lon, vals, 
        vmin=None, vmax=None, 
        clim=None,
        colorbar_label=None, 
        title=None,
):

    rad = np.deg2rad(90) - np.abs(lat)
    theta = lon 

    # contour data over the map.
    im = ax.contourf(theta, rad, vals, vmin=vmin, vmax=vmax)  
    cb = plt.colorbar(im, ax=ax)
    ax.set_title(title) if title is not None else None


def wrap(dictarr):
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


def calc_efield(phi, lat, lon):
    """ Calculate E-field from potential phi, using partial derivative approach 
    Assume phi is in Volts, output will be Volts/metre
    """
    
    shapeE = (phi.shape[0] - 2, phi.shape[1] - 2)
    En = np.zeros(shapeE) 
    Ee = np.zeros(shapeE) 
    latE = np.zeros(shapeE) 
    lonE = np.zeros(shapeE) 
    wgs84 = nv.FrameE(name='WGS84')
    z = -780E3  # depth in m

    for lati in range(shapeE[0]):
        for loni in range(shapeE[1]):
            latE[lati, loni] = lat[lati+1, loni+1]
            lonE[lati, loni] = lon[lati+1, loni+1]

            pointA = wgs84.GeoPoint(
                latitude=lat[lati, loni+1], 
                longitude=lon[lati, loni+1], 
                z=z, degrees=True,
            )
            pointB = wgs84.GeoPoint(
                latitude=lat[lati+2, loni+1], 
                longitude=lon[lati+2, loni+1], 
                z=z, degrees=True,
            )
            pointC = wgs84.GeoPoint(
                latitude=lat[lati+1, loni], 
                longitude=lon[lati+1, loni], 
                z=z, degrees=True,
            )
            pointD = wgs84.GeoPoint(
                latitude=lat[lati+1, loni+2], 
                longitude=lon[lati+1, loni+2], 
                z=z, degrees=True,
            )
    
            gc_distLat_m = pointA.delta_to(pointB).length
            gc_distLon_m = pointC.delta_to(pointD).length
            
            En[lati, loni] = (phi[lati, loni + 1] - phi[lati + 2, loni + 1]) / gc_distLat_m  # direction of increasing AACGM lat
            Ee[lati, loni] = (phi[lati + 1, loni + 2] - phi[lati + 1, loni]) / gc_distLon_m  # direction of eastward AAGCM lon

    return En, Ee, latE, lonE
            

def grad_sph(sfld, theta, phi, rad=65E3):
    # take gradient of scalar field in theta/phi directions
    d_sfld = np.gradient(sfld, [theta, phi])
    d_sfld_dtheta = 1 / r * d_sfld[0]
    d_sfld_dphi = 1 / (r * np.sin(theta)) * d_sfld[1]
     

def get_data(fname, dati=None):
    dset = ncread_vars(fname)

    # if no date/time requested assume that we're doing the whole
    # AMPERE netcdf file
    # otherwise, do just the requested time
     
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


if __name__ == '__main__':
    main()
