#!/usr/bin/python3 
# Plot data for a single day  
# - Plots electric potential maps with ion drift velocity overlaid 
# - Uses MAGNETIC coordinates converted to polar (r,theta)
# - Loads DMSP data for the same day and computes the ion drift velocity 
#   and overlays a quiver plot on those computed from the potential maps  
# - Save the plot to a file 
# FIXME
# 1. Problems with the shape of the DMSP data; of order 13k events, 1D arrays; 
#    Not quite sure how this maps to the ion drifts derived for (r,th) ~ (5,20) 
#    size from the electric potential maps.  Maybe this means we are fundamentally 
#    wrong when making the quiver plot(s).   
# 2. DMSP data has intrinsically np.MaskedArrays in the dictionary, which makes 
#    applying cuts on them a bit cumbersome; these data were previously downloaded 
#    from cedar.openmadrigal.org.  Trying to use pysatMadrigal, the data is _always_ empty. 
#    Not sure if this is correlated to issue #1.  

import os
from datetime import datetime, timedelta
import pysat
import pyIGRF
import ppigrf
import scipy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import nvector as nv
import my_utils
import electromagnetic as em 

import graph

#_______________________________________________________________________________

if __name__=='__main__':

    debug     = False
    plot_dmsp = True 

    graph.setGraphParameters(fontSize=15)

    # settings 
    # format is year, month, day, hr, min, sec 
    time = datetime(2019,3,1,0,10,0)
    b_is_up         = True    # is B assumed to have only up component (in ENU coordinates)    
    z               = 350E+3  # height for calculations [m] 
    noon_lat_offset = 55.     # radius at which to plot the noon dot 
    N_div           =  10     # bin size for contour plot color bar axis  
    V_min           = -200    # electric potential limit for color bar axis [kV]  
    V_max           =  200    # electric potential limit for color bar axis [kV]   

    # load data
    fpath   = my_utils.get_file_path(time) 
    data    = my_utils.read_file(fpath)
    outpath = './plots/v_ion_{0}-{1}-{2}_{3}-{4}-{5}.png'.format(time.year,time.month,time.day,time.hour,time.minute,time.second)
    # get north magnetic pole location in geographic coordinates
    np_idx      = data['MLAT (AACGM)'] == 90
    np_latlon   = [ data['Geographic Latitude'][np_idx][0], data['Geographic Longitude'][np_idx][0] ]

    # compute the E field 
    E_dict = em.calculate_E_mix(data,z)
    En     = E_dict['En'][1:,:].flatten() 
    Ee     = E_dict['Ee'][1:,:].flatten() 
    # get the magnetic field 
    B0,Be,Bn,Bu = em.get_B(E_dict['glat'],E_dict['glon'],time,z,b_is_up)
    if debug: 
        print("mean Ee = {0:.3E} V/m".format(np.mean(Ee)))
        print("mean En = {0:.3E} V/m".format(np.mean(En)))
        print("mean B0 = {0:.3E} T".format(np.mean(B0)))  
        print("mean Be = {0:.3E} T".format(np.mean(Be)))  
        print("mean Bn = {0:.3E} T".format(np.mean(Bn)))  
        print("mean Bu = {0:.3E} T".format(np.mean(Bu)))  

    # compute ion drift velocity 
    v_ion_mag,v_ion_dir,v_ion_e,v_ion_n = em.get_exb_drift_velocity(E_dict,B0,Bu,E_dict['glat'],E_dict['glon'])

    if debug:
        print("mean v_ion = {0:.3E} m/s".format(np.mean(v_ion_mag))) 

    if plot_dmsp:
        # load DMSP data
        sat_id      = 15
        inpath_dmsp = my_utils.get_file_path(time,'DMSP',sat_id)
        dmsp        = my_utils.read_file(inpath_dmsp)
        dmsp        = my_utils.process_dmsp(dmsp,mlat_cut=60,np_latlon=np_latlon)
        # compute ion drift in (N,E,U) and prepare for quiver plot  
        # mlat_dmsp           = np.deg2rad(90) - np.deg2rad(dmsp['mlat'])
        # mlon_dmsp           = np.deg2rad(dmsp['mlong'])
        # radius_dmsp,theta_dmsp = my_utils.get_radius_and_theta(np.unique(mlat_dmsp),np.unique(mlon_dmsp)) 
        radius_dmsp         = np.deg2rad(90) - np.deg2rad(dmsp['mlat'][:-1]) 
        theta_dmsp          = np.deg2rad(dmsp['mlong'][:-1]) 
        vi_n_dmsp,vi_e_dmsp = my_utils.bearing_magnitude_to_North_East(np.deg2rad(dmsp['vi_dirn_MAG']),dmsp['vi_mag']) 
        print('vi_n_dmsp = {0}, vi_e_dmsp = {1}'.format(vi_n_dmsp.shape,vi_e_dmsp.shape))  
        x_dmsp,y_dmsp       = my_utils.pol2cart_vec(radius_dmsp,theta_dmsp,-vi_n_dmsp,vi_e_dmsp)
        print('x_dmsp    = {0}, y_dmsp    = {1}'.format(x_dmsp.shape,y_dmsp.shape))  

    # make a plot 
    fig,ax = plt.subplots(1,1,subplot_kw={'projection':'polar'})

    # get (magnetic) radius and theta 
    radius,theta = my_utils.get_radius_and_theta(np.unique(data['MLAT (AACGM)']),np.unique(data['MLON (AACGM)']))
    print('radius = {0}, theta = {1}'.format(radius.shape,theta.shape)) 

    levels = np.linspace(V_min,V_max,N_div) 
    im     = ax.contourf(theta,radius,data['Potential'][::-1,:-1],levels=levels) 
    yticks = 90 - np.rad2deg(ax.get_yticks()) 
    ax.set_yticklabels(['%.0f' % y for y in yticks])  
    cbar = plt.colorbar(im,ax=ax)
    cbar.set_label('Electric Potential [kV]') 

    # quiver ion drift plot [our values] 
    x,y = my_utils.pol2cart_vec(radius,theta,v_ion_n[::-1,:],-v_ion_e[::-1,:]) 
    # print('x = {0}, y = {1}'.format(x.shape,y.shape)) 
    th,rd = np.meshgrid(theta,radius) 
    # print('th = {0}, rd = {1}'.format(th.shape,rd.shape)) 
    ax.quiver(th,rd,x,y,color='black') 
    
    if plot_dmsp:
        # guessing we need a mesh grid for the plot? 
        # th_dmsp,rd_dmsp = np.meshgrid(theta_dmsp,radius_dmsp) 
        print('theta_dmsp = {0}, radius_dmsp = {1}, x_dmsp = {2}, y_dmsp = {3}'.format(theta_dmsp.shape,radius_dmsp.shape,x_dmsp.shape,y_dmsp.shape))  
        ax.quiver(theta_dmsp,radius_dmsp,x_dmsp,y_dmsp,color='m') 

    # local noon dot 
    noon_glon = my_utils.local_noon(time) 
    noon_mlon = my_utils.glon_to_mlon(data['Geographic Longitude'][-1,:],data['MLON (AACGM)'][-1,:],noon_glon) 
    noon_rad,noon_theta = my_utils.get_radius_and_theta(noon_lat_offset+0.5,noon_mlon) 
    
    ax.plot(noon_theta,noon_rad,'.r',markersize=20) 
    ax.set_rmax(np.deg2rad(90-noon_lat_offset))
    ax.set_title(time.ctime())  

    plt.show()
    # plt.savefig(outpath)
    # print('Figure saved to: {0}'.format(outpath)) 
    # plt.close()

