# utility functions 

import os
import errno
import netCDF4
import pprint
import pvlib 
import scipy
import numpy as np 
import nvector as nv 
import pysat
import pysatMadrigal as pysatMad
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
    print("[filter_data]: {0}/{1} events passed the cut {2} in range {3}".format(cntr,N,var_name,var))  
    return out 
#_______________________________________________________________________________
def process_dmsp(dmsp:dict,mlat_cut:float,np_latlon=None): 
    # some preprocessing steps for DMSP data
    dmsp        = filter_data(dmsp,'mlat',np.array([mlat_cut]))
    # FIXME: another way to cut on data. This doesn't work because of masked arrays... 
    # breakpoint() 
    # dmsp        = dmsp[(dmsp['mlat']>mlat_cut)&np.isfinite(dmsp['hor_ion_v'])] 
    # calculate bearings
    # note: dimension of output array bearings is 1 less than the inputs since
    # we're taking a difference of the input arrays 
    time        = np.array(dmsp['timestamps']) 
    lat,lon,alt = np.array(dmsp['gdlat']), np.array(dmsp['glon']), np.array(dmsp['gdalt']) 
    _,bearings    = calculate_bearings(lat[:-1],lon[:-1],alt[:-1],lat[1:],lon[1:])
    # calculate velocity direction
    dmsp['vi_mag'] = dmsp['hor_ion_v']
    vi_mag = np.array(dmsp['vi_mag']) 
    dmsp['vi_dirn_geo']  = calculate_velocity_direction(lat[:-1],lon[:-1],time[:-1],bearings,vi_mag[:-1]) 
    dmsp['vi_mag_model'] = np.ones(len(dmsp))*np.nan
    # adjust vi_mag length; since we'll use it with vi_dirn_MAG to get v_ion in (N,E,U)
    # they must have the same length 
    dmsp['vi_mag'] = vi_mag[:-1] 
    # conversion to MAG drift directions if north pole is specified  
    if np_latlon is not None: 
        _,np_bearing = calculate_bearings(lat,lon,alt,np.ones(lat.shape)*np_latlon[0],np.ones(lat.shape)*np_latlon[1])
        dmsp['vi_dirn_MAG'] = dmsp['vi_dirn_geo'] + np_bearing[:-1]
    return dmsp 
#_______________________________________________________________________________
def load_dmsp_2(sat,stime,etime,dl_dmsp=False,usr='AlexChartier',pw='alex.chartier@jhuapl.edu'):
    """ Load DMSP data using pysatMadrigal. Optionally download """
    # Create PySat instrument
    dmsp = pysat.Instrument(
        inst_module=pysatMad.instruments.dmsp_ivm,
        inst_id=sat,
        clean_level='clean',
    )
    if dl_dmsp:
        dmsp.download(stime,etime,user=usr,password=pw)

    dmsp.load(date=stime)
    if dmsp.data.empty:
        print('No data for %s on %s' % (sat, stime.strftime('%Y %b %d %H:%M')))
        return []

    return dmsp
#_______________________________________________________________________________
def process_dmsp_2(dmsp,dec_rate=1,mlat_cutoff=60,np_latlon=None):
    """ Process DMSP """
    assert mlat_cutoff > 0, 'SH filtering not implemented yet'
    try:
        dmsp_vals = dmsp.data[::dec_rate] # decimate
    except:
        return []
    dmsp_vals  = dmsp_vals[(dmsp_vals['mlat'] > mlat_cutoff) & np.isfinite(dmsp_vals['hor_ion_v'])]
    dmsp_times = dmsp_vals.index
    lats,lons,alts = dmsp_vals['gdlat'], dmsp_vals['glon'], dmsp_vals['gdalt']
    brng_deg = calculate_bearings(lats[:-1], lons[:-1], alts[:-1], lats[1:], lons[1:])
    dmsp_vals = dmsp_vals[:-1]
    # add new keys to the dictionary 
    dmsp_vals['vi_mag']        = dmsp_vals['hor_ion_v']
    dmsp_vals['vi_dirn_geo']   = calculate_velocity_direction(lats[:-1],lons[:-1],dmsp_times[:-1],brng_deg,dmsp_vals['vi_mag'])
    dmsp_vals['vi_mag_model']  = np.ones(len(dmsp_vals)) * np.nan
    # dmsp_vals['bearings_ecef'] = brng_ecef  
    dmsp_vals['bearings']  = brng_deg
    # dmsp_vals['grid_pt_ecef']  = grid_pt_ecef 
    # add conversion to MAG drift directions if north pole is provided
    if np_latlon:
        np_bearings = calculate_bearings(lats, lons, alts,
            np.ones(lats.shape) * np_latlon[0], np.ones(lats.shape) * np_latlon[1])
        dmsp_vals['vi_dirn_MAG'] = dmsp_vals['vi_dirn_geo'] + np_bearings[:-1]

    return dmsp_vals # also return the bearings 
#_______________________________________________________________________________
def bearing_magnitude_to_North_East(vi_dir:np.array,vi_mag:np.array):
    # convert vector magnitude and direction to (N,E) components
    ND = len(vi_dir)
    NM = len(vi_mag)  
    if(ND!=NM):
        msg = 'len(vi_dir) = {0}, len(vi_mag) = {1}'.format(ND,NM)
        raise ValueError(msg)
    vi_N = np.cos(vi_dir)*vi_mag 
    vi_E = np.sin(vi_dir)*vi_mag 
    return vi_N,vi_E 
#_______________________________________________________________________________
def calculate_bearings(lat:np.array,lon:np.array,alt:np.array,latbs:np.array,lonbs:np.array):
    # compute bearings using the reference ellipsoid [deg/km] 
    wgs84 = nv.FrameE(name='WGS84')
    assert len(lat)==len(latbs),'Assuming lat and latbs are equal length arrays'
    bearing_deg = np.zeros(len(lat))*np.nan
    # bearing     = np.zeros(len(lat))*np.nan  
    # grid_pt     = np.zeros(len(lat))*np.nan # this is the starting point of the vector for the bearings   
    for idx, latb in enumerate(latbs): 
        depth   = (-1.)*alt[idx]*1E+3 # convert to meters 
        point_B = wgs84.GeoPoint(latitude=lat[idx],longitude=lon[idx]  ,z=depth,degrees=True) 
        point_A = wgs84.GeoPoint(latitude=latb    ,longitude=lonbs[idx],z=depth,degrees=True) # TODO: create a vector of these in ECEF (and return)
        p_AB_N  = point_A.delta_to(point_B) # we want the bearing at point A                  # TODO: create a vector of these in ECEF (and return)
        # grid_pt[idx]     = point_A.to_ecef_vector() 
        # bearing[idx]     = p_AB_N.to_ecef_vector()  
        bearing_deg[idx] = p_AB_N.azimuth_deg
    return bearing_deg 
#_______________________________________________________________________________
def get_bearings_ecef(lat:np.array,lon:np.array,alt:np.array,yaw=0,pitch=0,roll=0,deg=True): 
    # convert geodetic (lat,lon,z) to ECEF (x,y,z) using nvector
    wgs84 = nv.FrameE(name='WGS84')
    z     = -1.*alt
    N     = len(lat) - 1
    XYZ   = np.zeros((N,3))
    UVW   = np.zeros((N,3))    
    UVWR  = np.zeros((N,3))    
    for i in range(N):
        p  = wgs84.GeoPoint(latitude=lat[i]  ,longitude=lon[i]  ,z=z[i]  ,degrees=deg)
        p2 = wgs84.GeoPoint(latitude=lat[i+1],longitude=lon[i+1],z=z[i+1],degrees=deg)
        # compute the relative distance to the next (lat,lon,z) point
        # note: can effectively turn this into a bearing angle by calling .azimuth_deg
        brng      = p.delta_to(p2)
        # rotate the vector
        frame_B   = nv.FrameB(brng.to_nvector(),yaw=yaw,pitch=pitch,roll=roll,degrees=deg)
        p_brng_B  = frame_B.Pvector(np.r_[-1,0,0].reshape((-1,1))) # argument is a unit vector along the x axis (along axis of vehicle)
        # convert to ECEF frame
        p_ecef        = p.to_ecef_vector()
        brng_ecef     = brng.to_ecef_vector()
        brng_ecef_rot = p_brng_B.to_ecef_vector()
        # fill the output
        XYZ[i,:]  = p_ecef.pvector[:].flatten()
        UVW[i,:]  = brng_ecef.pvector[:].flatten()
        UVWR[i,:] = brng_ecef_rot.pvector[:].flatten()

    return XYZ,UVW
#_______________________________________________________________________________
def calculate_velocity_direction(lat:np.array,lon:np.array,time:np.array,bearing:np.array,vel:np.array,debug=False): 
    # DMSP horizontal ion drifts provided as follows: 
    # hor_ion_v = horizontal ion velocity (positive = sunward) [m/s] 
    # Here, we determine the direction
    if debug:
        print('lat: {0}'.format(lat.shape))
        print('lon: {0}'.format(lon.shape))
        print('time: {0}'.format(time.shape))
        print('bearing: {0}'.format(bearing.shape))
        print('vel: {0}'.format(vel.shape))

    vel_dir       = np.zeros(vel.shape)
    ephem_df      = pvlib.solarposition.get_solarposition(time,lat,lon)
    solaz         = ephem_df['azimuth']# [:-1] 
    bearing_solaz = zero_360(bearing-solaz)
    bearing_90    = zero_360(bearing+90)  
    bearing_270   = zero_360(bearing-90) 

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
def bearing_to_xy_comp(bearing:np.array): 
    # convert bearing angle to (x,y) angular components 
    NB = len(bearing) 
    x  = np.zeros(NB) 
    y  = np.zeros(NB)
    th = np.zeros(NB)
    for i in range(len(bearing)):
        bb = bearing[i]
        if bb<0:
            bb += 360. 
        # find vector angle 
        if bb>0 and bb<=90:
            # NE quadrant 
            theta = 90. - bb 
        elif bb>90 and bb<=180:
            # SE quadrant  
            theta = 360 - (bearing[i]-90.)
        elif bb>180 and bb<=270: 
            # SW quadrant
            theta = 270. - bb 
        elif bb>270 and bb<=360: 
            # NW quadrant   
            theta = 180 - (270.-bb) 
        x[i]  = np.cos(theta)
        y[i]  = np.sin(theta) 
        th[i] = theta  
    return x,y,th
#_______________________________________________________________________________
def geodetic_latlonz_to_ecef(lat:np.array,lon:np.array,z:np.array,deg=True,scale='m'):
    # convert geodetic (lat,lon,z) to ECEF (x,y,z) using nvector
    sf = 1.
    if scale=='km': 
        sf = 1e-3 
    wgs84 = nv.FrameE(name='WGS84')
    N = len(lat) - 1
    XYZ = np.zeros((N, 3))
    UVW = np.zeros((N, 3))
    for i in range(N):
        p  = wgs84.GeoPoint(latitude=lat[i]  ,longitude=lon[i]  ,z=z[i]  ,degrees=deg)
        p2 = wgs84.GeoPoint(latitude=lat[i+1],longitude=lon[i+1],z=z[i+1],degrees=deg)
        # compute the relative distance to the next (lat,lon,z) point
        brng      = p.delta_to(p2)     # can effectively turn this into a bearing angle by calling .azimuth_deg
        # convert to ECEF frame
        p_ecef    = p.to_ecef_vector()
        brng_ecef = brng.to_ecef_vector()
        # fill the output
        XYZ[i, :] = sf*p_ecef.pvector[:].flatten()      # starting point of vector (x0,y0,z0)  
        UVW[i, :] = sf*brng_ecef.pvector[:].flatten()   # magnitude of components/direction of vector (u,v,w) 
    return XYZ,UVW
# #_______________________________________________________________________________
# def geodetic_latlonz_to_ecef(lat:np.array,lon:np.array,z:np.array,deg=True):
#     # convert geodetic (lat,lon,z) to ECEF (x,y,z) using nvector  
#     wgs84  = nv.FrameE(name='WGS84') 
#     N      = len(lat)
#     x_ecef = np.zeros(N) 
#     y_ecef = np.zeros(N) 
#     z_ecef = np.zeros(N) 
#     for i in range(N): 
#         p         = wgs84.GeoPoint(latitude=lat[i],longitude=lon[i],z=z[i],degrees=deg) 
#         p_ecef    = p.to_ecef_vector()
#         pv        = p_ecef.pvector.ravel().tolist()  
#         x_ecef[i] = pv[0]  
#         y_ecef[i] = pv[1]  
#         z_ecef[i] = pv[2] 
#     return x_ecef,y_ecef,z_ecef 
#_______________________________________________________________________________
def drift_velocity_enu2ecef(lat_deg:np.array,lon_deg:np.array,vi_dirn:np.array,vi_mag:np.array,sf=1.0):
    # convert v_ion from ENU to ECEF coordinates 
    N1 = len(lat_deg) 
    N2 = len(vi_dirn) 
    if(N1!=N2):
        msg = 'ERROR! len(lat_deg) = {0}, len(vi_dirn) = {1}'.format(N1,N2) 
        raise ValueError(msg)
    # get east and north components of the velocity
    vi_dirn_rad = np.deg2rad(vi_dirn)  
    vi_E = vi_mag*np.sin(vi_dirn_rad) 
    vi_N = vi_mag*np.cos(vi_dirn_rad)
    uvw = np.zeros((N1-1,3))
    for i in range(0,N1-1):
        # rotation matrix to go from ENU -> ECEF 
        Lambda = np.deg2rad(lon_deg[i])  
        Phi    = np.deg2rad(lat_deg[i])  
        M      = np.matrix([[-np.sin(Lambda),-np.sin(Phi)*np.cos(Lambda),np.cos(Phi)*np.cos(Lambda)],
                            [np.cos(Lambda) ,-np.sin(Phi)*np.sin(Lambda),np.cos(Phi)*np.sin(Lambda)],
                            [0,np.cos(Phi),np.sin(Lambda)]])
        # vector in ENU
        q_enu  = np.matrix([[vi_E[i]],[vi_N[i]],[0]])  
        # conversion to ECEF 
        q_ecef = M.dot(q_enu)
        # test output 
        # print('lat = {0:.7f}, lon = {1:.7f}, q(enu) = {2}, q(ecef) = {3}'.format(lat_deg[i],lon_deg[i],np.asarray(q_enu).ravel(),np.asarray(q_ecef).ravel()))
        # fill output  
        uvw[i,:] = sf*np.asarray(q_ecef).ravel()  
    return uvw  
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
    '''for matplotlib polar quiver'''
    # FIXME: Not implemented properly...
    # NDR = len(dr)
    # NDT = len(dt)  
    # NTH = len(theta)  
    # if(NDR!=NDT or NTH!=NDR or NTH!=NDT):
    #     msg = 'dr = {0}, dt = {1}, theta = {2}'.format(dr.shape,dt.shape,theta.shape)
    #     raise ValueError(msg)
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
def ecef_to_pvector(x:np.array):
    # convert ecef to p vector
    N = len(x)  
    v = np.zeros((N,3)) 
    for i in range(0,N):
        v[i,:] = x[i].pvector[:].flatten() 
    return v 
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
        # prefix    = './data/pymix_{0:04d}{1}{2:02d}'.format(time.year,month_str,time.day)
        prefix    = './data/pymix/{0}{1}'.format(month_str.lower(),time.year-2000)
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
