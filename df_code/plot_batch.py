#!/usr/bin/python3 
# Plot data in discreet time steps between a start and stop time  
# - Plots electric potential maps with ion drift velocity overlaid 
# - Uses MAGNETIC coordinates converted to polar (r,theta) 
# - Saves each plot to a file 

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

    debug = False

    graph.setGraphParameters(fontSize=15)

    # settings 
    # format is year, month, day, hr, min, sec 
    start_time = datetime(2019,3,1,0,10,0)
    end_time   = datetime(2019,3,2,0,0 ,0)   
    time_step  = timedelta(minutes=2) 
    b_is_up         = True    # is B assumed to have only up component (in ENU coordinates)    
    z               = 350E+3  # height for calculations [m] 
    noon_lat_offset = 55.     # radius at which to plot the noon dot 
    N_div           =  10     # bin size for contour plot color bar axis  
    V_min           = -200    # electric potential limit for color bar axis [kV]  
    V_max           =  200    # electric potential limit for color bar axis [kV]   

    time = start_time 
    while(time<end_time): 
        # load data
        fpath   = my_utils.get_file_path(time) 
        data    = my_utils.read_file(fpath)
        outpath = './plots/v_ion_{0}-{1}-{2}_{3}-{4}-{5}.png'.format(time.year,time.month,time.day,time.hour,time.minute,time.second)
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

        # make a plot 
        fig,ax = plt.subplots(1,1,subplot_kw={'projection':'polar'})

        # get (magnetic) radius and theta 
        radius,theta = my_utils.get_radius_and_theta(np.unique(data['MLAT (AACGM)']),np.unique(data['MLON (AACGM)']))

        levels = np.linspace(V_min,V_max,N_div) 
        im     = ax.contourf(theta,radius,data['Potential'][::-1,:-1],levels=levels) 
        yticks = 90 - np.rad2deg(ax.get_yticks()) 
        ax.set_yticklabels(['%.0f' % y for y in yticks])  
        cbar = plt.colorbar(im,ax=ax)
        cbar.set_label('Electric Potential [kV]') 

        # quiver ion drift plot 
        x,y = my_utils.pol2cart_vec(radius,theta,v_ion_n[::-1,:],-v_ion_e[::-1,:]) 
        ax.quiver(theta,radius,x,y) 

        # local noon dot 
        noon_glon = my_utils.local_noon(time) 
        noon_mlon = my_utils.glon_to_mlon(data['Geographic Longitude'][-1,:],data['MLON (AACGM)'][-1,:],noon_glon) 
        noon_rad,noon_theta = my_utils.get_radius_and_theta(noon_lat_offset+0.5,noon_mlon) 
        
        ax.plot(noon_theta,noon_rad,'.r',markersize=20) 
        ax.set_rmax(np.deg2rad(90-noon_lat_offset))
        ax.set_title(time.ctime())  

        plt.savefig(outpath)
        print('Figure saved to: {0}'.format(outpath)) 
        plt.close()
       
        # set up for next time step 
        time += time_step  


