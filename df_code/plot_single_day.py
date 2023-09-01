#!/usr/bin/python3 
# Plot data for a single day  
# - Plots electric potential maps with ion drift velocity overlaid 
# - Uses MAGNETIC coordinates converted to polar (r,theta)
# - Loads DMSP data for the same day and computes the ion drift velocity 
#   and overlays a quiver plot on those computed from the potential maps  
# - Save the plot to a file 
# FIXME
# 1. Show that the bearings are perpendicular to the velocities measured by the satellite
#    a. Plot the satellite bearing angles on the quiver plot.  Likely need to convert bearing 
#       angles into equivalent magnetic (lat,lon) angles...   
#    b. Come up with a 'dot product' type evaluation to show they are perpendicular?  
# 2. Plot the satellite and our ion drift velocities in ECEF coordinates for an easier comparison 
#    a. Use nvector for the transformation  
# 3. ['done'] Quiver plot for DMSP and our data [focus on timestamp 2019-03-02T11-16-00Z]  
#    a. [DONE] Compute a test path for what the satellite path should look like (c.f., figure 4 in paper).
#       Extract (r,th) and make a quiver plot. Is this close to what we have in our data?  
#    b. [DONE] Divide DMSP data up into small time steps (~10 min) and make the quiver plot. 
#       Choose the time step based on getting a full path over the northern hemisphere.    
#    c. [DONE] Compute the model ion drift velo at the DMSP trajectory; use median time of that 
#       trajectory when evaluating our model.  Better to stay in the DMSP coordinates.   

import os
from datetime import datetime, timedelta
import aacgmv2
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
    dec_factor      = 2 # decimation factor for DMSP data
    mlat_cut        = 60 
    width           = 0.002 

    graph.setGraphParameters(fontSize=15)

    # settings 
    # format is year, month, day, hr, min, sec 
    # time = datetime(2019,3,1,0,10,0)
    time = datetime(2019,3,2,11,16,0)
    b_is_up         = True     # is B assumed to have only up component (in ENU coordinates)    
    z0              = -780E+3  # height for calculations [m] 
    noon_lat_offset = 55.      # radius at which to plot the noon dot 
    N_div           =  10      # bin size for contour plot color bar axis  
    V_min           = -90      # electric potential limit for color bar axis [kV]  
    V_max           =  45      # electric potential limit for color bar axis [kV]   

    # load data
    fpath   = my_utils.get_file_path(time) 
    data    = my_utils.read_file(fpath)
    outpath = './plots/v_ion_{0}-{1}-{2}_{3}-{4}-{5}.png'.format(time.year,time.month,time.day,time.hour,time.minute,time.second)
    # get north magnetic pole location in geographic coordinates
    np_idx      = data['MLAT (AACGM)'] == 90
    np_latlon   = [ data['Geographic Latitude'][np_idx][0], data['Geographic Longitude'][np_idx][0] ]

    # make a plot 
    fig,ax = plt.subplots(1,1,subplot_kw={'projection':'polar'})

    # get (magnetic) radius and theta 
    radius,theta = my_utils.get_radius_and_theta(np.unique(data['MLAT (AACGM)']),np.unique(data['MLON (AACGM)']))
    print('radius = {0}, theta = {1}'.format(radius.shape,theta.shape)) 

    # plot the electric potential 
    levels = np.linspace(V_min,V_max,N_div) 
    im     = ax.contourf(theta,radius,data['Potential'][::-1,:-1],levels=levels) 
    yticks = 90 - np.rad2deg(ax.get_yticks()) 
    ax.set_yticklabels(['%.0f' % y for y in yticks])  
    cbar = plt.colorbar(im,ax=ax)
    cbar.set_label('Electric Potential [kV]') 
    
    fig2 = plt.figure() 
    ax2  = fig2.add_subplot(projection='3d') 
    fig3 = plt.figure() 
    ax3  = fig3.add_subplot(projection='3d') 
    # fig4,ax4 = plt.subplots(1,1)
    # fig5 = plt.figure()
    # ax5 = fig5.add_subplot(projection='polar') 

    if plot_dmsp:
        # load DMSP data
        dmsp = {} 
        sat_str = ''
        for sat in sats:
            print("--------------------")  
            print("Loading DMSP data for satellite {0}...".format(sat)) 
            dmsp_data = my_utils.load_dmsp_2(sat,time,time+timedelta(days=1),dl_dmsp)
            print("Processing...") 
            dmsp[sat] = my_utils.process_dmsp_2(dmsp_data,dec_factor,mlat_cut,np_latlon)
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
                vi_N,vi_E = my_utils.bearing_magnitude_to_North_East(np.deg2rad(dmsp[sat]['vi_dirn_MAG'][tidx]),dmsp[sat]['vi_mag'][tidx],)
                x,y  = my_utils.pol2cart_vec(rad,theta,-vi_N,vi_E)
                ax.quiver(theta,rad,x,y,color='orange',width=width)
                print("--> Done!")
                print('Computing ion drift velocity at DMSP (glat, glon) for sat {0}'.format(sat))
                # now compute v_ion using DMSP (glat,glon)
                # print(dmsp[sat].keys()) 
                # print(dmsp[sat]) 
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
                bearings    = dmsp[sat]['bearings'][tidx].to_frame().loc[:,'bearings'].to_numpy()  
                vi_dirn_geo = dmsp[sat]['vi_dirn_geo'][tidx].to_frame().loc[:,'vi_dirn_geo'].to_numpy()
                # use average time 
                tt_ns = np.zeros(len(yy))
                for i in range(len(yy)):
                    tt_ns[i] = datetime(yy[i],mm[i],dd[i],hh[i],MIN[i],sec[i]).timestamp()
                tt_mean = np.mean(tt_ns)
                time_mean = datetime.fromtimestamp(tt_mean) 
                # plot the bearings
                xb,yb,thb   = my_utils.bearing_to_xy_comp(bearings)
                rb = np.sqrt( np.power(xb,2.) + np.power(yb,2.) ) 
                # convert bearings to magnetic coordinates
                mblat  = np.zeros(len(bearings)) 
                mbrngs = np.zeros(len(bearings)) 
                for i in range(len(bearings)):
                    mlat,mlon,malt = aacgmv2.get_aacgm_coord(glat[i],bearings[i],alt[i]/1E+3,time_mean)
                    mblat[i]  = mlat 
                    mbrngs[i] = mlon
                ax.plot(np.deg2rad(bearings),rad,color='black' ,markersize=6)      # plot the bearings in GEOGRAPHIC coordinates
                ax.plot(np.deg2rad(mbrngs)  ,np.deg2rad(mblat),color='magenta',markersize=6) # plot the bearings in MAGNETIC coordinates
                # for i in range(len(bearings)):
                #     print("bearing = {0:.2f} deg, theta = {1:.2f} deg, r = {2:.2f}, x = {3:.2f}, y = {4:.2f}".format(bearings[i],thb[i],rb[i],xb[i],yb[i]))  
                # make a plot of bearing angles vs vi_dirn_geo (should be perp.) 
                # ax4.plot(bearings,vi_dirn_geo)  
                z0 = (-1)*np.mean(alt) # -780 km is about the height of the satellite (- => nvector expects depth; -depth => above surface)  
                print('tt:   {0}'.format(tt_ns.shape)) 
                print('glat: {0}'.format(glat.shape))
                print('glon: {0}'.format(glon.shape))
                print('alt:  {0}'.format(alt.shape) )
                print('bearings:  {0}'.format(bearings.shape) )
                print('mean time: {0}'.format(time_mean))  
                E_dict = em.calculate_E_mix(data,z0)
                En     = E_dict['En'][1:,:].flatten() 
                Ee     = E_dict['Ee'][1:,:].flatten() 
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
                x,y   = my_utils.pol2cart_vec(rad,theta,v_ion_n,v_ion_e)
                # ax.quiver(theta,rad,x,y,color='blue',width=width)
                ax2.scatter(glat,glon,v_ion_n,label='ExB' )
                ax2.scatter(glat,glon,vi_N   ,label='DMSP')
                ax3.scatter(glat,glon,v_ion_e,label='ExB' )
                ax3.scatter(glat,glon,vi_E   ,label='DMSP')
                print("--> Done!") 

    # local noon dot 
    noon_glon = my_utils.local_noon(time) 
    noon_mlon = my_utils.glon_to_mlon(data['Geographic Longitude'][-1,:],data['MLON (AACGM)'][-1,:],noon_glon) 
    noon_rad,noon_theta = my_utils.get_radius_and_theta(noon_lat_offset+0.5,noon_mlon) 
    
    ax.plot(noon_theta,noon_rad,'.r',markersize=20) 
    ax.set_rmax(np.deg2rad(90-noon_lat_offset))
 
    title = '{0} {1}'.format(sat_str,time.ctime())
    ax.set_title(title)  

    fig2.suptitle('North Component')
    ax2.set_xlabel('GLAT [deg]')
    ax2.set_ylabel('GLON [deg]')
    ax2.set_zlabel(r'$v_{ion}^{n}$ [m/s]')
    ax2.legend(loc='best')

    fig3.suptitle('East Component')
    ax3.set_xlabel('GLAT [deg]')
    ax3.set_ylabel('GLON [deg]')
    ax3.set_zlabel(r'$v_{ion}^{e}$ [m/s]')
    ax3.legend(loc='best')

    # ax4.set_title('Bearings vs vi_dirn_geo') 
    # ax4.set_xlabel('Bearing Angle [deg]') 
    # ax4.set_ylabel('vi_dirn_geo [deg]') 
    # ax4.tick_params(axis='both')
    # ax4.xaxis.grid(True,which='both',linestyle='--')
    # ax4.yaxis.grid(True,which='both',linestyle='--')

    plt.show()
    # plt.savefig(outpath)
    # print('Figure saved to: {0}'.format(outpath)) 
    # plt.close()

