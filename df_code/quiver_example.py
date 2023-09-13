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

if __name__=='__main__':

    graph.setGraphParameters(fontSize=17,lineWidth=1) 

    plt.rcParams["figure.figsize"]    = [10,8]
    plt.rcParams["figure.autolayout"] = True

    # x = np.arange(-1,1,0.1) 
    # y = np.arange(-1,1,0.1) 

    # X,Y = np.meshgrid(x,y) 
    # z   = my_func(X,Y) 
    # dX,dY = np.gradient(z)  
    # 
    # f = plt.figure()
    # 
    # ax = f.add_subplot(polar=False)
    # ax.quiver(X,Y,dX,dY)
    # ax.set_aspect('equal') 
    # ax.set_xlabel('x [arb]')  
    # ax.set_ylabel('y [arb]')  

    # defining a rotation matrix about the z axis (in ECEF)  
    th = np.deg2rad(90.)  
    Mz = np.matrix([[np.cos(th),-np.sin(th),0],
                    [np.sin(th),np.cos(th),0],
                    [0,0,1]
                   ])
    My = np.matrix([[np.cos(th),0,np.sin(th)],
                    [0,1,0],
                    [-np.sin(th),0,np.cos(th)]
                   ])

 
    # MM = nv.rotation.R2xyz(M)
 
    lat = np.arange(60, 75)
    lon = np.zeros(lat.shape)
    alt = np.ones(lat.shape) * 800E3
    z = -alt
    deg = True
 
    # convert geodetic (lat,lon,z) to ECEF (x,y,z) using nvector
    wgs84 = nv.FrameE(name='WGS84')
    N    = len(lat) - 1
    XYZ  = np.zeros((N, 3))
    UVW  = np.zeros((N, 3))
    XYZR = np.zeros((N, 3))
    UVWR = np.zeros((N, 3))
    for i in range(N):
        p  = wgs84.GeoPoint(latitude=lat[i],longitude=lon[i],z=z[i],degrees=deg)
        p2 = wgs84.GeoPoint(latitude=lat[i+1],longitude=lon[i+1],z=z[i+1],degrees=deg)
        # compute the relative distance to the next (lat,lon,z) point  
        brng      = p.delta_to(p2)     # note: can effectively turn this into a bearing angle by calling .azimuth_deg
        brng_cpy  = brng     
        # convert to ECEF frame 
        p_ecef    = p.to_ecef_vector()
        brng_ecef = brng.to_ecef_vector()
        # rotate the vector 
        # temp       = p_ecef.pvector[:].flatten() 
        # temp_r     = M.dot(temp) 
        # p_ecef_rot = nv.FrameE(name='WGS84').ECEFvector(temp_r)  
        start     = p_ecef.pvector[:].flatten()
        temp      = brng_ecef.pvector[:].flatten()
        temp_r    = Mz.dot(temp)
        # temp_r[0,1] -= temp_r[0,1]  
        # temp_r    = My.dot(temp_r.reshape(3,1))
        # temp_r    = np.ravel(temp_r)  
        # # temp_r    = temp_r[0,:].flatten()
        # print(type(temp))
        # print(type(temp_r))
        print('start: {0}, end: {1}, end after rotation: {2}'.format(start,temp,temp_r))  
        brng_ecef_rot = nv.FrameE(name='WGS84').ECEFvector(temp_r) 
        # brng_ecef_rot = brng_ecef.rotation.R_EL2n_E(M) # doesn't work  
        # fill the output  
        XYZ[i,:]  = p_ecef.pvector[:].flatten()
        UVW[i,:]  = brng_ecef.pvector[:].flatten()
        # XYZR[i,:] = p_ecef_rot.pvector[:].flatten()
        UVWR[i,:] = brng_ecef_rot.pvector[:].flatten()
 
    # breakpoint()
    ax = plt.figure().add_subplot(projection='3d')
    ax.quiver(*XYZ.T,*UVW.T)
    ax.quiver(*XYZ.T,*UVWR.T,color='red')

    ax.set_xlabel('x') 
    ax.set_ylabel('y') 
    ax.set_zlabel('z') 
    ax.set_xlim([0,4E+6])
    ax.set_ylim([-4,4])
    ax.set_zlim([0,7E+6])
 
    plt.show()
