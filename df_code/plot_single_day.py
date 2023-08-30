#!/usr/bin/python3 
# Plot data for a single day  
# - Plots electric potential maps with ion drift velocity overlaid 
# - Uses MAGNETIC coordinates converted to polar (r,theta)
# - Loads DMSP data for the same day and computes the ion drift velocity 
#   and overlays a quiver plot on those computed from the potential maps  
# - Save the plot to a file 
# FIXME
# 1. Can make a quiver plot; but seeing effectively multiple trajectories 
#    for the DMSP data.  Not consistent with expectations. 
#    Next Steps: [focus on timestamp 2019-03-02T11-16-00Z]  
#    a. [DONE] Compute a test path for what the satellite path should look like (c.f., figure 4 in paper).
#       Extract (r,th) and make a quiver plot. Is this close to what we have in our data?  
#    b. [DONE] Divide DMSP data up into small time steps (~10 min) and make the quiver plot. 
#       Choose the time step based on getting a full path over the northern hemisphere.    
#    c. Compute the model ion drift velo at the DMSP trajectory; use median time of that 
#       trajectory when evaluating our model.  Better to stay in the DMSP coordinates.   
# 2. DMSP data has intrinsically np.MaskedArrays in the dictionary, which makes 
#    applying cuts on them a bit cumbersome; these data were previously downloaded 
#    from cedar.openmadrigal.org.  Trying to use pysatMadrigal, the data is _always_ empty. 
#    Not sure if this is correlated to issue #1. 

import os
from datetime import datetime, timedelta
import pysat
import pysatMadrigal as pysatMad
import pyIGRF
import ppigrf
import scipy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import nvector as nv

import my_utils
import graph
import electromagnetic as em 

#_______________________________________________________________________________

if __name__=='__main__':

    debug           = False
    plot_map_quiver = False # plot quiver using E-field map coordinates 
    plot_dmsp       = True
    dl_dmsp         = False 
    sats            = ['f16','f17','f18']  # satellites (DMSP) 
    dec_rate        = 1 # decimation rate for DMSP data
    mlat_cut        = 60 
    width           = 0.002 

    graph.setGraphParameters(fontSize=15)

    # settings 
    # format is year, month, day, hr, min, sec 
    # time = datetime(2019,3,1,0,10,0)
    time = datetime(2019,3,2,11,16,0)
    b_is_up         = True    # is B assumed to have only up component (in ENU coordinates)    
    z               = 350E+3  # height for calculations [m] 
    noon_lat_offset = 55.     # radius at which to plot the noon dot 
    N_div           =  10     # bin size for contour plot color bar axis  
    V_min           = -90    # electric potential limit for color bar axis [kV]  
    V_max           =  45    # electric potential limit for color bar axis [kV]   

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
    if plot_map_quiver:
        x,y = my_utils.pol2cart_vec(radius,theta,v_ion_n[::-1,:],-v_ion_e[::-1,:]) 
        th,rd = np.meshgrid(theta,radius) 
        print('x = {0}, y = {1}'.format(x.shape,y.shape)) 
        print('th = {0}, rd = {1}'.format(th.shape,rd.shape)) 
        ax.quiver(th,rd,x,y,color='blue',width=width) 

    sat_str = ''
    
    fig2 = plt.figure() 
    ax2  = fig2.add_subplot(projection='3d') 
    fig3 = plt.figure() 
    ax3  = fig3.add_subplot(projection='3d') 

    if plot_dmsp:
        # load DMSP data
        dmsp = {} 
        for sat in sats:
            print("--------------------")  
            print("Loading DMSP data for satellite {0}...".format(sat)) 
            dmsp_data = my_utils.load_dmsp_2(sat,time,time+timedelta(days=1),dl_dmsp)
            print("Processing...") 
            dmsp[sat] = my_utils.process_dmsp_2(dmsp_data,dec_rate,mlat_cut,np_latlon)
            # get valid time indices for each satellite 
            tidx = np.abs(dmsp[sat].index - time) < timedelta(minutes=10) / 2
            if np.sum(tidx) == 0:
                print('No data in the time interval {0} to {1}.  Moving on...'.format(time,time+timedelta(days=1)))
                continue
            else: 
                sat_str += '{0} '.format(sat) 
                # create quiver plot 
                print("Plotting...") 
                rad   = np.deg2rad(90) - np.deg2rad(dmsp[sat]['mlat'][tidx])
                theta = np.deg2rad(dmsp[sat]['mlong'][tidx])
                vi_N, vi_E = my_utils.bearing_magnitude_to_North_East(np.deg2rad(dmsp[sat]['vi_dirn_MAG'][tidx]), dmsp[sat]['vi_mag'][tidx],)
                x, y  = my_utils.pol2cart_vec(rad,theta,-vi_N,vi_E)
                ax.quiver(theta,rad,x,y,color='m',width=width)
                print("--> Done!")
                print('Computing ion drift velocity at DMSP (glat, glon) for sat {0}'.format(sat))
                # now compute v_ion using DMSP (glat,glon)
                # print(dmsp[sat].keys()) 
                # DMSP position; have to extract it from the pandas dataframe and convert to np arrays 
                yy   = dmsp[sat]['year'][tidx].to_frame().loc[:,'year'].to_numpy() 
                mm   = dmsp[sat]['month'][tidx].to_frame().loc[:,'month'].to_numpy() 
                dd   = dmsp[sat]['day'][tidx].to_frame().loc[:,'day'].to_numpy() 
                hh   = dmsp[sat]['hour'][tidx].to_frame().loc[:,'hour'].to_numpy() 
                MIN  = dmsp[sat]['min'][tidx].to_frame().loc[:,'min'].to_numpy() 
                sec  = dmsp[sat]['sec'][tidx].to_frame().loc[:,'sec'].to_numpy() 
                glat = dmsp[sat]['gdlat'][tidx].to_frame().loc[:,'gdlat'].to_numpy() 
                glon = dmsp[sat]['glon'][tidx].to_frame().loc[:,'glon'].to_numpy() 
                alt  = dmsp[sat]['gdalt'][tidx].to_frame().loc[:,'gdalt'].to_numpy()*1E+3 # convert to meters 
                # FIXME: use average time? 
                tt_ns = np.zeros(len(yy))
                for i in range(len(yy)):
                    tt_ns[i] = datetime(yy[i],mm[i],dd[i],hh[i],MIN[i],sec[i]).timestamp()
                tt_mean = np.mean(tt_ns)
                time_mean = datetime.fromtimestamp(tt_mean)  
                print('tt:   {0}'.format(tt_ns.shape)) 
                print('glat: {0}'.format(glat.shape))
                print('glon: {0}'.format(glon.shape))
                print('alt:  {0}'.format(alt.shape) )
                print('mean time: {0}'.format(time_mean))  
                # compute magnetic field  
                B0,Be,Bn,Bu = em.get_B(glat,glon,time_mean,alt,b_is_up)
                # compute ion drift velocity 
                v_ion_mag,v_ion_dir,v_ion_e,v_ion_n = em.get_exb_drift_velocity(E_dict,B0,Bu,glat,glon,0.1)
                # make quiver plot
                # for i in range(len(yy)):
                #     print("glat = {0:.3f}, glon = {1:.3f}, vi_N: dmsp = {2:.3f} ours = {3:.3f}, vi_E: dmsp = {4:.3f} ours = {5:.3f}".format(glat[i],glon[i],vi_N[i],v_ion_n[i],vi_E[i],v_ion_e[i])) 
                print("Plotting...") 
                mlat  = dmsp[sat]['mlat'][tidx].to_frame().loc[:,'mlat'].to_numpy() 
                mlon  = dmsp[sat]['mlong'][tidx].to_frame().loc[:,'mlong'].to_numpy() 
                rad   = np.deg2rad(90) - np.deg2rad(mlat)
                theta = np.deg2rad(mlon)
                x,y   = my_utils.pol2cart_vec(rad,theta,v_ion_n,-v_ion_e)
                ax.quiver(theta,rad,x,y,color='black',width=width)
                ax2.scatter(glat,glon,v_ion_n)
                ax2.scatter(glat,glon,vi_N)
                ax3.scatter(glat,glon,v_ion_e)
                ax3.scatter(glat,glon,vi_E)
                print("--> Done!") 
 
        # inpath_dmsp = my_utils.get_file_path(time,'DMSP',sat_id)
        # dmsp        = my_utils.read_file(inpath_dmsp)
        # dmsp        = my_utils.process_dmsp(dmsp,mlat_cut=60,np_latlon=np_latlon)
        # # compute ion drift in (N,E,U) and prepare for quiver plot  
        # radius_dmsp,theta_dmsp = my_utils.get_radius_and_theta(dmsp['mlat'][:-1],dmsp['mlong'][:-1]) 
        # vi_n_dmsp,vi_e_dmsp = my_utils.bearing_magnitude_to_North_East(np.deg2rad(dmsp['vi_dirn_MAG']),dmsp['vi_mag']) 
        # print('vi_n_dmsp = {0}, vi_e_dmsp = {1}'.format(vi_n_dmsp.shape,vi_e_dmsp.shape))  
        # x_dmsp,y_dmsp       = my_utils.pol2cart_vec(radius_dmsp,theta_dmsp,-vi_n_dmsp,vi_e_dmsp)
        # print('x_dmsp    = {0}, y_dmsp    = {1}'.format(x_dmsp.shape,y_dmsp.shape))  
    
    # if plot_dmsp:
    #     # guessing we need a mesh grid for the plot? 
    #     # th_dmsp,rd_dmsp = np.meshgrid(theta_dmsp,radius_dmsp) 
    #     print('theta_dmsp = {0}, radius_dmsp = {1}, x_dmsp = {2}, y_dmsp = {3}'.format(theta_dmsp.shape,radius_dmsp.shape,x_dmsp.shape,y_dmsp.shape))  
    #     ax.quiver(theta_dmsp,radius_dmsp,x_dmsp,y_dmsp,color='m') 
    #     # dummy curve
    #     # my_lat = np.array([61. ,70. ,80. ,90. ,80.,70. ,61, ])
    #     # my_lon = np.array([315.,300.,270.,240.,210,180.,150.])
    #     # my_vin = np.array([1.  ,1.  ,1.  ,1.  ,1. ,1.  ,1.  ])
    #     # my_vie = np.array([1.  ,1.  ,1.  ,1.  ,1. ,1.  ,1.  ])
    #     # my_rad,my_th = my_utils.get_radius_and_theta(my_lat,my_lon)
    #     # my_x,my_y    = my_utils.pol2cart_vec(my_rad,my_th,my_vin,my_vie) 
    #     # ax.quiver(my_th,my_rad,my_x,my_y,color='tab:orange') 

    # local noon dot 
    noon_glon = my_utils.local_noon(time) 
    noon_mlon = my_utils.glon_to_mlon(data['Geographic Longitude'][-1,:],data['MLON (AACGM)'][-1,:],noon_glon) 
    noon_rad,noon_theta = my_utils.get_radius_and_theta(noon_lat_offset+0.5,noon_mlon) 
    
    ax.plot(noon_theta,noon_rad,'.r',markersize=20) 
    ax.set_rmax(np.deg2rad(90-noon_lat_offset))
 
    title = '{0} {1}'.format(sat_str,time.ctime())
    ax.set_title(title)  

    ax2.set_xlabel('GLAT [deg]')
    ax2.set_ylabel('GLON [deg]')
    ax2.set_zlabel(r'$v_{ion}^{n}$ [m/s]')

    ax3.set_xlabel('GLAT [deg]')
    ax3.set_ylabel('GLON [deg]')
    ax3.set_zlabel(r'$v_{ion}^{e}$ [m/s]')

    plt.show()
    # plt.savefig(outpath)
    # print('Figure saved to: {0}'.format(outpath)) 
    # plt.close()

