#!/usr/bin/python3 
# Quiver plot example
# Rotates bearings through an angle of 270 deg about the vehicle's vertical axis 
# and plots the initial and rotated vectors.  

import numpy as np
import matplotlib.pyplot as plt
import nvector as nv

import graph

#_______________________________________________________________________________
if __name__=='__main__':

    graph.setGraphParameters(fontSize=17,lineWidth=1) 

    plt.rcParams["figure.figsize"]    = [10,8]
    plt.rcParams["figure.autolayout"] = True
 
    lat = np.arange(60, 75)
    lon = np.zeros(lat.shape)
    alt = np.ones(lat.shape) * 800E3
    z   = -alt
    deg = True
 
    # convert geodetic (lat,lon,z) to ECEF (x,y,z) using nvector
    wgs84 = nv.FrameE(name='WGS84')
    N     = len(lat) - 1
    XYZ   = np.zeros((N, 3))
    UVW   = np.zeros((N, 3))
    UVWR  = np.zeros((N, 3))
    for i in range(N):
        p  = wgs84.GeoPoint(latitude=lat[i],longitude=lon[i],z=z[i],degrees=deg)
        p2 = wgs84.GeoPoint(latitude=lat[i+1],longitude=lon[i+1],z=z[i+1],degrees=deg)
        # compute the relative distance to the next (lat,lon,z) point 
        # note: can effectively turn this into a bearing angle by calling .azimuth_deg
        brng      = p.delta_to(p2)     
        # rotate the vector
        frame_B   = nv.FrameB(brng.to_nvector(),yaw=270,pitch=0,roll=0,degrees=deg)  
        p_brng_B  = frame_B.Pvector(np.r_[-1,0,0].reshape((-1,1))) # argument is a unit vector along the x axis (along axis of vehicle)  
        # convert to ECEF frame 
        p_ecef        = p.to_ecef_vector()
        brng_ecef     = brng.to_ecef_vector()
        brng_ecef_rot = p_brng_B.to_ecef_vector()  
        # fill the output  
        XYZ[i,:]  = p_ecef.pvector[:].flatten()
        UVW[i,:]  = brng_ecef.pvector[:].flatten()
        UVWR[i,:] = brng_ecef_rot.pvector[:].flatten()
 
    # breakpoint()
    ax = plt.figure().add_subplot(projection='3d')
    ax.quiver(*XYZ.T,*UVW.T)
    ax.quiver(*XYZ.T,*UVWR.T,color='red')

    ax.set_xlabel('x') 
    ax.set_ylabel('y') 
    ax.set_zlabel('z') 
    # ax.set_xlim([0,4E+6])
    ax.set_ylim([-4,4])
    # ax.set_zlim([0,7E+6])
 
    plt.show()
