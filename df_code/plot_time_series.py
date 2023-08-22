#!/usr/bin/python3 
# Calculate the ion drift velocity and plot the N and E components vs time  
# TODO 
# 1. Properly account for out-of-range values for glat, glon 
# 2. Properly account for varying satellite height; currently using average height 
# 3. Extract up component of v_ion (?) DMSP has only up and horiz components.  
# 4. Apply a cut to accept data for |mlat| >= 60 deg 

import os
import sys
from datetime import datetime, timedelta
import pysat
import pyIGRF
import ppigrf
import scipy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates 
import nvector as nv
import my_utils
import electromagnetic as em 

import math_util, graph 
# from utildf import algorithm, constants, graph, math_util, units

#_______________________________________________________________________________

if __name__=='__main__':

    debug = True 

    graph.setGraphParameters(fontSize=16)

    # settings 
    # format is year, month, day, hr, min, sec 
    start_time      = datetime(2019,3,1,0,10,0)
    end_time        = datetime(2019,3,1,2,10,0)   
    # end_time        = datetime(2019,3,2,0,0 ,0)   
    time_step       = timedelta(minutes=2) 
    b_is_up         = True    # is B assumed to have only up component (in ENU coordinates)    
    z               = 350E+3  # height for calculations [m] 
    # plot limits 
    y_min           = -500 
    y_max           =  500 

    # get north magnetic pole location in geographic coordinates 
    fpath       = my_utils.get_file_path(start_time)
    data        = my_utils.read_file(fpath)
    np_idx      = data['MLAT (AACGM)'] == 90 
    np_latlon   = [ data['Geographic Latitude'][np_idx][0], data['Geographic Longitude'][np_idx][0] ]
    
    # load DMSP data
    sat_id      = 15
    time_dmsp   = datetime(2019,3,1,0,0,0)  
    inpath_dmsp = my_utils.get_file_path(time_dmsp,'DMSP',sat_id)
    dmsp        = my_utils.read_file(inpath_dmsp)
    dmsp        = my_utils.process_dmsp(dmsp,mlat_cut=60,np_latlon=np_latlon) 

    # prepare for plotting 

    t_dmsp    = np.array( dmsp['timestamps'] ) 
    glat_dmsp = np.array( dmsp['gdlat']      )
    glon_dmsp = np.array( dmsp['glon']       )
    v_u_dmsp  = np.array( dmsp['vert_ion_v'] ) # FIXME: this is up!  
    v_e_dmsp  = np.array( dmsp['hor_ion_v']  )
    z_dmsp    = np.array( dmsp['gdalt']      )*1E+3 # convert from km to m  

    # remove NaNs 
    good_idx = np.argwhere(~np.isnan(v_u_dmsp)) 
    # print(good_idx) 
    t_dmsp    = np.take(t_dmsp,good_idx)  
    v_u_dmsp  = np.take(v_u_dmsp,good_idx)  
    v_e_dmsp  = np.take(v_e_dmsp,good_idx) 
    glat_dmsp = np.take(glat_dmsp,good_idx)  
    glon_dmsp = np.take(glon_dmsp,good_idx) 
    z_dmsp    = np.take(z_dmsp,good_idx)  

    v_mag_dmsp = np.sqrt( v_u_dmsp**2 + v_e_dmsp**2 ) 

    # convert timestamp to a datetime array
    dt = [] 
    for i in range(0,len(t_dmsp)): 
        dt.append(datetime.fromtimestamp(int(t_dmsp[i]))) 
    t_dmsp = np.array(dt)  

    print('--- DMSP Data ---') 
    print(t_dmsp.shape)
    print(glat_dmsp.shape)  
    print(glon_dmsp.shape)  
    print(v_u_dmsp.shape) 
    print(v_e_dmsp.shape) 
    print(z_dmsp.shape) 

    # also get ECEF versions
    v_x_dmsp = np.zeros(len(t_dmsp))  
    v_y_dmsp = np.zeros(len(t_dmsp))  
    v_z_dmsp = np.zeros(len(t_dmsp))  
    for i in range(0,len(t_dmsp)):
        v_e    = np.squeeze(v_e_dmsp)[i] 
        v_u    = np.squeeze(v_u_dmsp)[i]
        lat    = np.squeeze(glat_dmsp)[i] 
        lon    = np.squeeze(glon_dmsp)[i] 
        # print("{0},{1},{2}".format(t_dmsp[i],v_e,v_u))
        q_enu  = np.matrix([[v_e],[0],[v_u]])
        q_ecef = my_utils.enu2ecef(lat,lon,q_enu)
        v_x_dmsp[i] = q_ecef[0,0] 
        v_y_dmsp[i] = q_ecef[1,0] 
        v_z_dmsp[i] = q_ecef[2,0] 
 
    # output lists for plotting 
    t_out     = [] 
    v_mag_out = [] 
    v_n_out   = [] 
    v_e_out   = []
    v_x_out   = [] 
    v_y_out   = [] 
    v_z_out   = [] 
    glat_out  = [] 
    glon_out  = []
    z_out     = []  
 
    time = start_time 
    while(time<end_time): 
        # load data
        fpath   = my_utils.get_file_path(time)
        print("Processing file: {0}".format(fpath)) 
        data    = my_utils.read_file(fpath)
        # compute the E field 
        # note: the arrays in the output dictionary have the general shape (north,east)
        E_dict = em.calculate_E_mix(data,z)
        En     = E_dict['En'][1:,:].flatten() 
        Ee     = E_dict['Ee'][1:,:].flatten() 
        # get the DMSP satellite position based on timestamp  
        # ts     = datetime.timestamp(time)
        glat   = math_util.linearInterp(time,t_dmsp,glat_dmsp) 
        glon   = math_util.linearInterp(time,t_dmsp,glon_dmsp)
        z      = math_util.linearInterp(time,t_dmsp,z_dmsp   )  
        # get the magnetic field
        B0,Be,Bn,Bu = em.get_B(glat,glon,time,z,b_is_up)
        if debug: 
            print("    satellite {0} glat = {1:.5f} deg, glon = {2:.5f} deg, z = {3:.5f} km".format(sat_id,glat[0],glon[0],z[0]/1E+3)) 
            print("    mean Ee = {0:.3E} V/m".format(np.mean(Ee)))
            print("    mean En = {0:.3E} V/m".format(np.mean(En)))
            print("    mean B0 = {0:.3E} T".format(np.mean(B0)))  
            print("    mean Be = {0:.3E} T".format(np.mean(Be)))  
            print("    mean Bn = {0:.3E} T".format(np.mean(Bn)))  
            print("    mean Bu = {0:.3E} T".format(np.mean(Bu)))  
        # compute ion drift velocity 
        v_ion_mag,v_ion_dir,v_ion_e,v_ion_n = em.get_exb_drift_velocity(E_dict,B0,Bu,glat,glon,1)
        # convert from ENU to ECEF
        lat = np.squeeze(glat) 
        lon = np.squeeze(glon) 
        q_enu  = np.matrix([[v_ion_e],[v_ion_n],[0]])
        q_ecef = my_utils.enu2ecef(lat,glon,q_enu) 
        v_ion_x = np.squeeze(q_ecef[0,0])   
        v_ion_y = np.squeeze(q_ecef[1,0])   
        v_ion_z = np.squeeze(q_ecef[2,0]) 
        v_ecef_mag = np.sqrt( v_ion_x**2 + v_ion_y**2 + v_ion_z**2 )
        if debug:
            print("    {0},{1},{2}".format(v_ion_x,v_ion_y,v_ion_z))
            print("    mean |v_ion| = {0:.3E} m/s ({1:.3E})".format(np.mean(v_ion_mag),v_ecef_mag)) 
        # build the output arrays 
        t_out.append(time)
        v_mag_out.append(v_ion_mag) 
        v_n_out.append(v_ion_n)  
        v_e_out.append(v_ion_e)
        v_x_out.append(v_ion_x) 
        v_y_out.append(v_ion_y) 
        v_z_out.append(v_ion_z) 
        glat_out.append(glat[0])  
        glon_out.append(glon[0]) 
        z_out.append(z[0])  
        # move to next time 
        time += time_step  

    # convert to numpy arrays 
    t_out     = np.array(t_out) 
    v_mag_out = np.array(v_mag_out) 
    v_n_out   = np.array(v_n_out)  
    v_e_out   = np.array(v_e_out)  
    v_x_out   = np.array(v_x_out)  
    v_y_out   = np.array(v_y_out)  
    v_z_out   = np.array(v_z_out)  
    glat_out  = np.array(glat_out)  
    glon_out  = np.array(glon_out)  
    z_out     = np.array(z_out)  

    title_str    = 'Ion Drift Velocity [glat, glon change with satellite {0} location]'.format(sat_id)
    x_axis_title = 'Time' 
    y_axis_title = 'Drift Velocity [m/s]'
    date_format  = '%Y-%m-%d\n%H:%M:%S' 

    # build the figures
    NROW = 2
    NCOL = 1
 
    # north  
    fig_n,ax_n = plt.subplots(nrows=NROW,ncols=NCOL,sharex=False)
    fig_n.set_figheight(10)
    fig_n.set_figwidth(14)

    ax_n[0].plot(t_dmsp,v_u_dmsp,label=r'DMSP $v_{ion,up}$') 
    ax_n[0].plot(t_out ,v_n_out ,label=r'$v_{ion,n}$')
    # labels 
    ax_n[0].set_title('Ion Drift Velocity [North]')  
    ax_n[0].set_ylabel(y_axis_title)
    ax_n[1].set_ylabel('Percent Difference [%]')
    ax_n[1].set_xlabel(x_axis_title)

    for i in range(0,NROW): 
        ax_n[i].legend(loc='best') 
        # ax_n[i].set_xlabel(x_axis_title)
        # ax_n[i].set_ylim([y_min,y_max])
        ax_n[i].xaxis.set_major_formatter(mdates.DateFormatter(date_format)) 
        ax_n[i].tick_params(axis='both')  
        ax_n[i].xaxis.grid(True,which='both',linestyle='--') 
        ax_n[i].yaxis.grid(True,which='both',linestyle='--')
    # fig_n.autofmt_xdate() # rotate the date so it looks nice 

    # east 
    fig_e,ax_e = plt.subplots(nrows=NROW,ncols=NCOL,sharex=False)
    fig_e.set_figheight(10)
    fig_e.set_figwidth(14)

    ax_e[0].plot(t_dmsp,v_e_dmsp,label=r'DMSP $v_{ion,horiz}$') 
    ax_e[0].plot(t_out ,v_e_out ,label=r'$v_{ion,e}$') 
    # labels
    ax_e[0].set_title('Ion Drift Velocity [East]') 
    ax_e[0].set_ylabel(y_axis_title)
    ax_e[1].set_ylabel('Percent Difference [%]')
    ax_e[1].set_xlabel(x_axis_title)

    for i in range(0,NROW): 
        ax_e[i].legend(loc='best') 
        ax_e[i].set_xlabel(x_axis_title)
        # ax_n[i].set_ylim([y_min,y_max])
        ax_e[i].xaxis.set_major_formatter(mdates.DateFormatter(date_format)) 
        ax_e[i].tick_params(axis='both')  
        ax_e[i].xaxis.grid(True,which='both',linestyle='--') 
        ax_e[i].yaxis.grid(True,which='both',linestyle='--')
    # fig_e.autofmt_xdate() # rotate the date so it looks nice 

    # magnitude
    fig_m,ax_m = plt.subplots(nrows=NROW,ncols=NCOL,sharex=False)
    fig_m.set_figheight(10)
    fig_m.set_figwidth(14)

    ax_m[0].plot(t_dmsp,v_mag_dmsp,label=r'DMSP $\sqrt{v_{up}^{2} + v_{horiz}^{2}}$') 
    ax_m[0].plot(t_out ,v_mag_out ,label=r'$\sqrt{v_{n}^{2} + v_{e}^{2}}$')
    # labels 
    ax_m[0].set_title('Ion Drift Velocity Magnitude') 
    ax_m[0].set_ylabel(y_axis_title)
    ax_m[1].set_ylabel('Percent Difference [%]')
    ax_m[1].set_xlabel(x_axis_title)

    for i in range(0,NROW): 
        ax_m[i].legend(loc='best') 
        # ax_m[i].set_xlabel(x_axis_title)
        # ax_m[i].set_ylim([y_min,y_max])
        ax_m[i].xaxis.set_major_formatter(mdates.DateFormatter(date_format)) 
        ax_m[i].tick_params(axis='both')  
        ax_m[i].xaxis.grid(True,which='both',linestyle='--') 
        ax_m[i].yaxis.grid(True,which='both',linestyle='--')
    # fig_m.autofmt_xdate() # rotate the date so it looks nice 

    # ECEF 
    # x 
    fig_x,ax_x = plt.subplots(nrows=NROW,ncols=NCOL,sharex=False)
    fig_x.set_figheight(10)
    fig_x.set_figwidth(14)

    ax_x[0].plot(t_dmsp,v_x_dmsp,label=r'DMSP $v_{x}$') 
    ax_x[0].plot(t_out ,v_x_out ,label=r'$v_{x}$')
    # labels 
    ax_x[0].set_title('Ion Drift Velocity ECEF x Component') 
    ax_x[0].set_ylabel(y_axis_title)
    ax_x[1].set_ylabel('Percent Difference [%]')
    ax_x[1].set_xlabel(x_axis_title)

    for i in range(0,NROW): 
        ax_x[i].legend(loc='best') 
        ax_x[i].xaxis.set_major_formatter(mdates.DateFormatter(date_format)) 
        ax_x[i].tick_params(axis='both')  
        ax_x[i].xaxis.grid(True,which='both',linestyle='--') 
        ax_x[i].yaxis.grid(True,which='both',linestyle='--')

    # y 
    fig_y,ax_y = plt.subplots(nrows=NROW,ncols=NCOL,sharex=False)
    fig_y.set_figheight(10)
    fig_y.set_figwidth(14)

    ax_y[0].plot(t_dmsp,v_y_dmsp,label=r'DMSP $v_{y}$') 
    ax_y[0].plot(t_out ,v_y_out ,label=r'$v_{y}$')
    # labels 
    ax_y[0].set_title('Ion Drift Velocity ECEF y Component') 
    ax_y[0].set_ylabel(y_axis_title)
    ax_y[1].set_ylabel('Percent Difference [%]')
    ax_y[1].set_xlabel(x_axis_title)

    for i in range(0,NROW): 
        ax_y[i].legend(loc='best') 
        ax_y[i].xaxis.set_major_formatter(mdates.DateFormatter(date_format)) 
        ax_y[i].tick_params(axis='both')  
        ax_y[i].xaxis.grid(True,which='both',linestyle='--') 
        ax_y[i].yaxis.grid(True,which='both',linestyle='--')

    # z 
    fig_z,ax_z = plt.subplots(nrows=NROW,ncols=NCOL,sharex=False)
    fig_z.set_figheight(10)
    fig_z.set_figwidth(14)

    ax_z[0].plot(t_dmsp,v_z_dmsp,label=r'DMSP $v_{z}$') 
    ax_z[0].plot(t_out ,v_z_out ,label=r'$v_{z}$')
    # labels 
    ax_z[0].set_title('Ion Drift Velocity ECEF z Component') 
    ax_z[0].set_ylabel(y_axis_title)
    ax_z[1].set_ylabel('Percent Difference [%]')
    ax_z[1].set_xlabel(x_axis_title)

    for i in range(0,NROW): 
        ax_z[i].legend(loc='best') 
        ax_z[i].xaxis.set_major_formatter(mdates.DateFormatter(date_format)) 
        ax_z[i].tick_params(axis='both')  
        ax_z[i].xaxis.grid(True,which='both',linestyle='--') 
        ax_z[i].yaxis.grid(True,which='both',linestyle='--')

    # diagnostics 
   
    # plot the satellite glat, glon vs time 
    fig_d,ax_d = plt.subplots(nrows=3,ncols=1,sharex=False)
    fig_d.set_figheight(10)
    fig_d.set_figwidth(16)
   
    ax_d[0].plot(t_dmsp,glat_dmsp,label='DMSP')  
    ax_d[0].plot(t_out ,glat_out ,label='Interpolated')  

    ax_d[1].plot(t_dmsp,glon_dmsp,label='DMSP') 
    ax_d[1].plot(t_out ,glon_out ,label='Interpolated') 

    ax_d[2].plot(t_dmsp,z_dmsp/1E+3,label='DMSP') 
    ax_d[2].plot(t_out ,z_out/1E+3 ,label='Interpolated') 
    # labels
    ax_d[0].set_title('Lattitude') 
    ax_d[0].set_ylabel('Lattitude [deg]')
    ax_d[1].set_title('Longitude') 
    ax_d[1].set_ylabel('Longitude [deg]')
    ax_d[2].set_title('Altitude') 
    ax_d[2].set_ylabel('Altitude [km]')

    for i in range(0,3): 
        # for j in range(0,2): 
        ax_d[i].legend(loc='best') 
        # ax_d[i].set_xlabel(x_axis_title)
        ax_d[i].xaxis.set_major_formatter(mdates.DateFormatter(date_format)) 
        ax_d[i].tick_params(axis='both')  
        ax_d[i].xaxis.grid(True,which='both',linestyle='--') 
        ax_d[i].yaxis.grid(True,which='both',linestyle='--')
    # fig_d.autofmt_xdate() # rotate the date so it looks nice 

    plt.subplots_adjust(hspace=0.5)  
    plt.show()  
