#!/usr/bin/python3
# Quiver plot example

import numpy as np
import matplotlib.pyplot as plt
import nvector as nv
import os 
import nc_utils

#import graph


def plot_quiver(lat, lon, alt, deg=True):
    """ plots vectors between adjacent points, and vectors perpendicular to that """
    z = -alt

    # convert geodetic (lat,lon,z) to ECEF (x,y,z) using nvector
    wgs84 = nv.FrameE(name='WGS84')
    N = len(lat) - 1
    XYZ = np.zeros((N, 3))
    UVW = np.zeros((N, 3))
    XYZR = np.zeros((N, 3))
    UVWR = np.zeros((N, 3))


    for i in range(N):
        p = wgs84.GeoPoint(
            latitude=lat[i], longitude=lon[i], z=z[i], degrees=deg)
        p2 = wgs84.GeoPoint(
            latitude=lat[i+1], longitude=lon[i+1], z=z[i+1], degrees=deg)

        # compute the relative distance to the next (lat,lon,z) point
        # note: can effectively turn this into a bearing angle by calling .azimuth_deg
        brng = p.delta_to(p2)

        # calculate the rotation
        if np.sum(lat) > 0:
            frame_B = nv.FrameB(brng.to_nvector(), yaw=90, pitch=0, roll=0, degrees=deg)
        else:
            frame_B = nv.FrameB(brng.to_nvector(), yaw=270, pitch=0, roll=0, degrees=deg)
        p_BC_B = frame_B.Pvector(np.r_[1E5, 0, 0].reshape((-1, 1)))

        # convert to ECEF frame
        p_ecef = p.to_ecef_vector()
        brng_ecef = brng.to_ecef_vector()
        brng_ecef_rot = p_BC_B.to_ecef_vector()

        # fill the output
        XYZ[i, :] = p_ecef.pvector[:].flatten()
        UVW[i, :] = brng_ecef.pvector[:].flatten()
        # XYZR[i,:] = p_ecef_rot.pvector[:].flatten()
        UVWR[i, :] = brng_ecef_rot.pvector[:].flatten()

    """ do the plotting """
    plt.rcParams["figure.figsize"] = [8, 8]
    plt.rcParams["figure.autolayout"] = True
    ax = plt.figure().add_subplot(projection='3d')
    ax.quiver(*XYZ.T, *UVW.T, color='m')
    ax.quiver(*XYZ.T, *UVWR.T, color='red')

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    #ax.set_xlim([0, 4E+6])
    #ax.set_ylim([-4, 4])
    #ax.set_zlim([0, 7E+6])

    """ plot a sphere """ 
    # Define sphere parameters
    radius = 6371E3
    theta = np.linspace(0, 2.*np.pi, 100)
    phi = np.linspace(0, np.pi, 100)

    # Convert to Cartesian coordinates
    x = radius * np.outer(np.cos(theta), np.sin(phi))
    y = radius * np.outer(np.sin(theta), np.sin(phi))
    z = radius * np.outer(np.ones(np.size(theta)), np.cos(phi))
    ax.plot_surface(x, y, z, color='grey')

    plt.show()


def gen_dummy_vals():
    """ Create some input for the plotting """
    lat = np.arange(60, 75)
    lon = np.zeros(lat.shape)
    alt = np.ones(lat.shape) * 800E3
    return lat, lon, alt


def load_dmsp(dmsp_fn, lat_cutoff=45):
    """ load the DMSP data, filter and segment it """
    dmsp_fn = os.path.expanduser(dmsp_fn)
    dmsp = nc_utils.ncread_vars(dmsp_fn)

    # Decimate 
    for k, v in dmsp.items():
        dmsp[k] = v[::10]

    # Filter in lat 
    goodidx = np.abs(dmsp['gdlat']) > lat_cutoff
    for k, v in dmsp.items():
        dmsp[k] = v[goodidx]
    
    # figure out the segments based on a diff of > 10 degs longitude
    segs = np.where(np.diff(np.sign(dmsp['gdlat'])) != 0)
    
    dmsp_segs = []
    i0 = 0
    for i in segs[0]:
        i+=1
        tmp = {}
        for k, v in dmsp.items():
            tmp[k] = v[i0:i]

        goodi = np.isfinite(tmp['ion_v_sat_left'])
        for k, v in tmp.items():
            tmp[k] = v[goodi]
        dmsp_segs.append(tmp)
        i0 = i
        
    return dmsp_segs


if __name__ == '__main__':
    dmsp_fn = 'data/dmsp/dms_ut_20190302_15.002.nc'  # note: UT dallas files required
    dmsp = load_dmsp(dmsp_fn)
    
    lat, lon, alt = gen_dummy_vals()
    #breakpoint()
    plot_quiver(dmsp[0]['gdlat'], dmsp[0]['glon'], dmsp[0]['gdalt']*1E3)
    #plot_quiver(lat, lon, alt)


