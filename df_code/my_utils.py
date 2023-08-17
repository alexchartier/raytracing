# utility functions 

import os
import errno
import netCDF4
import scipy
import numpy as np 
import nvector as nv 
from datetime import datetime 

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
