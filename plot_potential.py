import numpy as np
import h5py
import os
import datetime as dt
import matplotlib.pyplot as plt
import nvector as nv
import pandas as pd
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
):

    """
    Validate potential maps using DMSP cross-track drift data
    """

    """ Get the North Magnetic Pole location in geographic coordinates """
    pot_fname = pot_fname_fmt % (stime.year, stime.month, stime.day, stime.hour, stime.minute, stime.second)
    pot = nc_utils.ncread_vars(pot_fname)

    np_idx = pot['MLAT (AACGM)'] == 90
    np_latlon = [pot['Geographic Latitude'][np_idx][0], pot['Geographic Longitude'][np_idx][0]]

    """ Run through and process """
    time = stime
    while time < etime:

        """ load poetntial at time """ 
        pot_fname = pot_fname_fmt % (time.year, time.month, time.day, time.hour, time.minute, time.second)
        pot = nc_utils.ncread_vars(pot_fname)

        """  plot """
        plot_fname = plot_fname_fmt % (time.year, time.month, time.day, time.hour, time.minute, time.second)
        fig, ax = plt.subplots(1, 1,subplot_kw={'projection': 'polar'}) 

        rad, theta = get_rad_theta(
            np.unique(pot['MLAT (AACGM)']), np.unique(pot['MLON (AACGM)']),
        )
        print(theta)

        """ contour data over the map. """
        latlim = 55.
   
        theta = np.concatenate([theta, [np.abs(theta[0]),]])
        poten = pot['Potential'][::-1, :]
        im = ax.contourf(theta, rad, poten)
        yticks = 90 - np.rad2deg(ax.get_yticks())
        ax.set_yticklabels(['%1.0f' % y for y in yticks])
        cb = plt.colorbar(im, ax=ax)
    
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
        print('***** Saved to %s' % plot_fname)
        plt.close()

        time += timestep 


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


def zero_360(vals):
    """ set values between 0 - 360 degrees """

    if isinstance(vals, float):
        vals = np.array(vals)
    vals[vals >= 360] -= 360
    vals[vals < 0] += 360
    return vals


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




















