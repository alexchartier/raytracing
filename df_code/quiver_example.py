#!/usr/bin/python3 
# Quiver plot example 

import numpy as np
import matplotlib.pyplot as plt
import nvector as nv

import graph

#_______________________________________________________________________________
def get_magnitude(x:np.array): 
    return np.sqrt( np.sum( np.power(x,2) ) ) 
#_______________________________________________________________________________
def my_func(x:np.ndarray,y:np.ndarray): 
    f = x*np.exp(-(x**2+y**2)) 
    return f
#_______________________________________________________________________________
def get_vector_field():
    x = np.arange(-1,1,0.1) 
    y = np.arange(-1,1,0.1)
    # create grid  
    X,Y = np.meshgrid(x,y)
    # evaluate function 
    z   = my_func(X,Y) 
    dX,dY = np.gradient(z)  
    return X,Y,dX,dY 
#_______________________________________________________________________________

if __name__=='__main__':

    graph.setGraphParameters(fontSize=17,lineWidth=1) 

    plt.rcParams["figure.figsize"]    = [10,8]
    plt.rcParams["figure.autolayout"] = True

    # X,Y,dX,dY = get_vector_field() 
    # f  = plt.figure()
    # ax = f.add_subplot(polar=False)
    # ax.quiver(X,Y,dX,dY)
    # ax.set_aspect('equal') 
    # ax.set_xlabel('x [arb]')  
    # ax.set_ylabel('y [arb]')  

    # defining a rotation matrix about the z axis (in ECEF)  
    th = np.deg2rad(-90.)  
    Mz = np.matrix([[np.cos(th),-np.sin(th),0],
                    [np.sin(th),np.cos(th),0],
                    [0,0,1]
                   ])
 
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
    UVWR2 = np.zeros((N, 3))
    for i in range(N):
        p  = wgs84.GeoPoint(latitude=lat[i],longitude=lon[i],z=z[i],degrees=deg)
        p2 = wgs84.GeoPoint(latitude=lat[i+1],longitude=lon[i+1],z=z[i+1],degrees=deg)
        # compute the relative distance to the next (lat,lon,z) point 
        # note: can effectively turn this into a bearing angle by calling .azimuth_deg
        brng      = p.delta_to(p2)     
        # rotate the vector
        frame_B   = nv.FrameB(brng.to_nvector(),yaw=90,pitch=0,roll=0,degrees=deg)  
        p_brng_B  = frame_B.Pvector(np.r_[1,0,0].reshape((-1,1))) # argument is a unit vector along the x axis 
        # wrong way
        brng_ecef = brng.to_ecef_vector() 
        temp      = brng_ecef.pvector[:].flatten()
        temp_r    = Mz.dot(temp)
        brng_ecef_rot_wrong = nv.FrameE(name='WGS84').ECEFvector(temp_r)     
        # convert to ECEF frame 
        p_ecef        = p.to_ecef_vector()
        brng_ecef     = brng.to_ecef_vector()
        brng_ecef_rot = p_brng_B.to_ecef_vector()  
        # fill the output  
        XYZ[i,:]  = p_ecef.pvector[:].flatten()
        UVW[i,:]  = brng_ecef.pvector[:].flatten()
        UVWR[i,:] = brng_ecef_rot.pvector[:].flatten()
        UVWR2[i,:] = brng_ecef_rot_wrong.pvector[:].flatten()
 
    # breakpoint()
    ax = plt.figure().add_subplot(projection='3d')
    ax.quiver(*XYZ.T,*UVW.T)
    ax.quiver(*XYZ.T,*UVWR.T,color='red')
    ax.quiver(*XYZ.T,*UVWR2.T,color='green')

    ax.set_xlabel('x') 
    ax.set_ylabel('y') 
    ax.set_zlabel('z') 
    # ax.set_xlim([0,4E+6])
    ax.set_ylim([-4,4])
    # ax.set_zlim([0,7E+6])
 
    plt.show()
