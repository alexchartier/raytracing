import numpy as np
import scipy.interpolate


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


def brng_mag_to_N_E(vi_dir, vi_mag):
    """ go from vector magnitude and direction to N, E components """
    vi_N = np.cos(vi_dir) * vi_mag
    vi_E = np.sin(vi_dir) * vi_mag

    return vi_N, vi_E


def zero_24(vals):
    """ set values between 0 - 24 hours """
    if isinstance(vals, float):
        vals = np.array(vals)
    vals[vals >= 24] -= 24
    vals[vals < 0] += 24
    return vals


def zero_360(vals):
    """ set values between 0 - 360 degrees """
    if isinstance(vals, float):
        vals = np.array(vals)
    vals[vals >= 360] -= 360
    vals[vals < 0] += 360
    return vals


def calc_lt(datetimes, lons):
    """ Function to calculate the local times from list of UTs and longitudes """
    UT = np.array([t.hour + t.minute/60 for t in datetimes])
    LT = UT + lons / 360 * 24 
    LT = zero_24(LT)

    return LT


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



