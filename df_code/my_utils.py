# utility functions 

import os
import errno
import netCDF4
import pvlib 
import scipy
import numpy as np 
import nvector as nv 
from datetime import datetime 

#_______________________________________________________________________________
def filter_data(data:dict,var_name:str,var:np.array):
    # filter a dictionary based on the variable var_name in the range var = [min,max]  
    # list of keys
    keys = [] 
    for entry in data: 
        keys.append(entry)
    # variable we cut on 
    values = data[var_name]
    # output dictionary 
    out = {}
    # cut on data 
    cntr = 0
    NC = len(var)  
    N  = len(values) 
    for i in range(0,N):
        passed = False
        # check the cut condition based on size of the cut array 
        if(NC==1):
            if(values[i]>var[0]):
                passed = True
        elif(NC==2):
            if(values[i]>var[0] and values[i]<var[1]):
                passed = True 
        if passed: 
            # passed the cut, populate output dictionary 
            cntr = cntr + 1 
            for k in keys:
                vv = data[k][i]
                # print("{0}: {1} (i = {2})".format(k,vv,i)) 
                if k in out:
                    out[k] += [vv] 
                else: 
                    out[k]  = [vv] 
    print("{0}/{1} events passed the cut {2} in range {3}".format(cntr,N,var_name,var))  
    return out 
#_______________________________________________________________________________
def process_dmsp(dmsp:dict,mlat_cut:float,np_latlon=None): 
    # some preprocessing steps for DMSP data
    dmsp        = filter_data(dmsp,'mlat',np.array([mlat_cut])) 
    # calculate bearings
    time        = np.array(dmsp['timestamps']) 
    lat,lon,alt = np.array(dmsp['gdlat']), np.array(dmsp['glon']), np.array(dmsp['gdalt']) 
    bearings    = calculate_bearings(lat[:-1],lon[:-1],alt[:-1],lat[1:],lon[1:])
    # calculate velocity direction
    dmsp['vi_mag'] = dmsp['hor_ion_v']
    vi_mag = np.array(dmsp['vi_mag']) 
    dmsp['vi_dirn_geo']  = calculate_velocity_direction(lat[:-1],lon[:-1],time[:-1],bearings,vi_mag[:-1]) # FIXME: why do we have to remove an element here?
    dmsp['vi_mag_model'] = np.ones(len(dmsp))*np.nan
    # conversion to MAG drift directions if north pole is specified  
    if np_latlon is not None: 
        np_bearing = calculate_bearings(lat,lon,alt,np.ones(lat.shape)*np_latlon[0],np.ones(lat.shape)*np_latlon[1])
        dmsp['vi_dirn_MAG'] = dmsp['vi_dirn_geo'] + np_bearing[:-1]
    return dmsp 
#_______________________________________________________________________________
def bearing_magnitude_to_North_East(vi_dir:np.array,vi_mag:np.array):
    # convert vector magnitude and direction to (N,E) components 
    vi_N = np.cos(vi_dir)*vi_mag 
    vi_E = np.sin(vi_dir)*vi_mag 
    return vi_N,vi_E 
#_______________________________________________________________________________
def calculate_bearings(lat:np.array,lon:np.array,alt:np.array,latbs:np.array,lonbs:np.array):
    # compute bearings using the reference ellipsoid [deg/km] 
    wgs84 = nv.FrameE(name='WGS84')
    assert len(lat)==len(latbs),'Assuming lat and latbs are equal length arrays'
    bearing_deg = np.zeros(len(lat))*np.nan 
    for idx, latb in enumerate(latbs): 
        depth   = (-1.)*alt[idx]*1E+3 # convert to meters 
        point_B = wgs84.GeoPoint(latitude=lat[idx],longitude=lon[idx]  ,z=depth,degrees=True) 
        point_A = wgs84.GeoPoint(latitude=latb    ,longitude=lonbs[idx],z=depth,degrees=True) 
        p_AB_N  = point_A.delta_to(point_B) # we want the bearing at point A 
        bearing_deg[idx] = p_AB_N.azimuth_deg 
    return bearing_deg 
#_______________________________________________________________________________
def calculate_velocity_direction(lat:np.array,lon:np.array,time:np.array,bearing:np.array,vel:np.array): 
    # DMSP horizontal ion drifts provided as follows: 
    # hor_ion_v = horizontal ion velocity (positive = sunward) [m/s] 
    # Here, we determine the direction
    # FIXME: the line where we get bearing_solaz ends up with 1 less element than all previous entries, 
    #        and hence can't carry on with the calculation.  
    print('time = {0}'.format(time.shape))  
    print('lat = {0}'.format(lat.shape))  
    print('lon = {0}'.format(lon.shape))  
    print('vel = {0}'.format(vel.shape))  
    vel_dir       = np.zeros(vel.shape)
    print('vel_dir = {0}'.format(vel_dir.shape))
    ephem_df      = pvlib.solarposition.get_solarposition(time,lat,lon)
    print('ephem_df = {0}'.format(ephem_df.shape))  
    solaz         = ephem_df['azimuth']# [:-1] 
    print('solaz = {0}'.format(solaz.shape)) 
    print('bearing = {0}'.format(bearing.shape))  
    bearing_solaz = zero_360(bearing-solaz)
    print('bearing_solaz = {0}'.format(bearing_solaz.shape))  
    bearing_90    = zero_360(bearing+90)  
    print('bearing_90 = {0}'.format(bearing_90.shape))  
    bearing_270   = zero_360(bearing-90) 
    print('bearing_270 = {0}'.format(bearing_270.shape))  

    # case A: sun on the right, positive velocity 
    id = (bearing_solaz<180) & (vel>=0)
    vel_dir[id] = bearing_90[id] 
    # case B: sun on the left, positive velocity 
    id = (bearing_solaz>=180) & (vel>=0) 
    vel_dir[id] = bearing_270[id] 
    # case C: sun on the right, negative velocity 
    id = (bearing_solaz<180) & (vel<0) 
    vel_dir[id] = bearing_90[id] 
    # case D: sun on the left, negative velocity 
    id = (bearing_solaz>=180) & (vel<0) 
    vel_dir[id] = bearing_270[id] 

    return vel_dir
#_______________________________________________________________________________
def enu2ecef(lat_deg:float,lon_deg:float,q_enu:np.array):
    # convert vector q in the East-North-Up system to that in the  
    # Earth-Centered Earth-Fixed (ECEF) system
    Lambda = np.deg2rad(lon_deg)  
    Phi    = np.deg2rad(lat_deg)  
    M      = np.matrix([[-np.sin(Lambda),-np.sin(Phi)*np.cos(Lambda),np.cos(Phi)*np.cos(Lambda)],
                        [np.cos(Lambda) ,-np.sin(Phi)*np.sin(Lambda),np.cos(Phi)*np.sin(Lambda)],
                        [0,np.cos(Phi),np.sin(Lambda)]]) 
    q_ecef = M.dot(q_enu) 
    return q_ecef
#_______________________________________________________________________________
def pol2cart_vec(radius:np.ndarray,theta:np.ndarray,dr:np.ndarray,dt:np.ndarray): 
    '''matplotlib polar quiver'''
    # FIXME: Not implemented properly... 
    x = dr*np.cos(theta) - dt*np.sin(theta) 
    y = dr*np.sin(theta) + dt*np.cos(theta) 
    return x,y
#_______________________________________________________________________________
def mag_to_E_N(vi_mag:np.ndarray,vi_dir:np.ndarray):
    '''convert v_ion magnitude and direction to east and north components'''
    vi_E = vi_mag*np.sin(vi_dir) 
    vi_N = vi_mag*np.cos(vi_dir)
    return vi_E,vi_N  
#_______________________________________________________________________________
def glon_to_mlon(glon,mlon,glon_i):
    '''interpolate geographic longitude to AACGM longitude'''
    interp_obj = scipy.interpolate.interp1d(zero_360(glon),zero_360(mlon),fill_value='extrapolate')
    noon_mlon  = interp_obj(zero_360(glon_i))
    return noon_mlon
#_______________________________________________________________________________
def zero_360(vals):
    '''set values between 0-360 deg'''
    if isinstance(vals,float):
        vals = np.array(vals)
    vals[vals>=360] -= 360
    vals[vals<0]    += 360
    return vals
#_______________________________________________________________________________
def local_noon(time:datetime):
    '''Calculate local noon longitude from UTC time'''
    noon_lon = zero_360(-(time.hour/24 + time.minute/1440)*360 + 180 )
    return noon_lon
#_______________________________________________________________________________
def get_radius_and_theta(lat:float,lon:float):
    '''Convert latitude and longitude to (r,theta)'''
    radius = np.deg2rad(90.) - np.abs(np.deg2rad(lat))
    theta  = np.deg2rad(lon)
    return radius,theta
#_______________________________________________________________________________
def wrap(d_input:dict):
    ''' prepare the data for processing '''
    d_output = {}
    for k,arr in d_input.items():
        # add opposite longitude ring above the pole
        half_lon    = int(arr.shape[1]/2)
        lon_idx     = np.concatenate([np.arange(half_lon,arr.shape[1]),np.arange(half_lon)])
        arr_2       = np.concatenate([arr[1,[lon_idx,]],arr],axis=0)
        # add duplicate ring below the lowest latitude
        arr_3       = np.concatenate([arr_2,arr_2[[-1,],:]],axis=0)
        # add another longitude wrap before the first one
        arr_4       = np.concatenate([arr_3[:,[-2,]],arr_3],axis=1)
        d_output[k] = arr_4

    return d_output
#_______________________________________________________________________________
def get_file_path(time:datetime,type_name:str='ampere',sat_id:int=0):
    month_list = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    month_str  = month_list[time.month-1]
    if type_name=='ampere': 
        prefix    = './data/pymix_{0:04d}{1}{2:02d}'.format(time.year,month_str,time.day)
        file_name = 'ampere_mix_{0:04d}-{1:02d}-{2:02d}T{3:02d}-{4:02d}-{5:02d}Z.nc'.format(time.year,time.month,time.day,time.hour,time.minute,time.second)
    elif type_name=='DMSP' or type_name=='dmsp': 
        prefix    = './data/dmsp' 
        file_name = 'dms_{0:04d}{1:02d}{2:02d}_{3}s1.001.nc'.format(time.year,time.month,time.day,sat_id) 
    else: 
        msg = 'Invalid data set: {0}'.format(type_name) 
        raise ValueError(msg) 
    fpath      = '{0}/{1}'.format(prefix,file_name)
    return fpath 
#_______________________________________________________________________________
def read_file(inpath:str,print_info=False):
    if isinstance(inpath,str): 
        fn = os.path.expanduser(inpath)
        if not os.path.isfile(fn): 
            raise FileNotFoundError(errno.ENOENT,os.strerror(errno.ENOENT),fn) 
        data = netCDF4.Dataset(fn,'r',format='NETCDF4')
    elif isinstance(inpath,netCDF4._netCDF4.Dataset): 
        data = inpath 
    else:
        raise ValueError('inpath needs to be a string or netCDF data set') 

    out = {} 

    for key in data.variables.keys(): 
        out[key] = data.variables[key][...] 
    data.close() 

    if(print_info):
        print("-----------------------------------------------")
        print("file: {0}".format(inpath))
        print(out)
        print("-----------------------------------------------")

    return out
