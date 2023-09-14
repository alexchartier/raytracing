# E x B ion drift calculations 

from datetime import datetime
import numpy as np 
import scipy 
import my_utils
import nvector as nv 
import ppigrf 

#_______________________________________________________________________________
def get_B(lat:np.array,lon:np.array,time:datetime,height:float,b_is_up:bool): 
    '''get the magnetic field using the IGRF model. Use geodetic coordinate  
    function call to igrf, which accounts for the ellipsoidal shape of the earth 
    input: 
    - lat     = geographic latitude [deg] 
    - lon     = geographic longitude [deg]
    - time    = time of day 
    - height  = height above sea level [m] 
    - b_is_up = B assumed to have only an upwards component [true or false]  
    output: 
    - B0      = Magnitude of the magnetic field [T] 
    - Be      = East component of the magnetic field [T]  
    - Bn      = North component of the magnetic field [T]  
    - Bu      = Up component of the magnetic field [T]  
    '''
    height_km = height/1E+3 
    Be,Bn,Bu = ppigrf.igrf(lon,lat,height_km,time) # igrf needs the height in km 
    # convert from nT to T  
    Be = Be*1e-9 
    Bn = Bn*1e-9 
    Bu = Bu*1e-9 
    # compute the magnitude
    if b_is_up:
       # assume B has only up component  
       B0 = np.squeeze( np.sqrt(Bu**2) )
    else: 
       B0 = np.squeeze( np.sqrt( Be**2 + Bn**2 + Bu**2 ) )
    # squeeze components for output # FIXME: Why not squeeze Bn, Be?
    Bu = np.squeeze(Bu) 
    return B0,Be,Bn,Bu
#_______________________________________________________________________________
def calculate_E_mix(data:dict,z:float):
    ''' calculate the E-field from the MIX potential.  Take the gradient in
    geo lat/lon coordinates.  We also need AACGM lat/lon for comparison with
    other data sets
    input:
    - data = AMPERE data file containing the electric potential, glat, glon,
             mlat, mlon data
    output:
    - mix  = dataframe containing the computed E-field (north and east components),
             glat, glon, mlat, mlon coordinates
    '''

    # set up output dataframe
    mix  = { }

    # prepare the data
    data_wrap = my_utils.wrap(data)

    # extract the data
    Phi  = data_wrap['Potential']*1E+3 # convert kV to V
    glat = data_wrap['Geographic Latitude']
    glon = data_wrap['Geographic Longitude']

    # compute the E field
    mix['En'],mix['Ee'],mix['glat'],mix['glon'] = calculate_E(Phi,glat,glon,z)

    # save the magnetic coordinates
    mix['mlat'] = data['MLAT (AACGM)'][:,:-1]
    mix['mlon'] = data['MLON (AACGM)'][:,:-1]

    return mix
#_______________________________________________________________________________
def calculate_E(Phi:np.ndarray,glat:np.ndarray,glon:np.ndarray,z:float):
    ''' calculate the E-field from the potential Phi.  Use the partial 
    derivative approach.
    input:
    - Phi    = Electric potential [V]
    - glat   = Geographic latitude [deg]
    - glon   = Geographic longitude [deg]
    output:
    - En     = Electric field in the north direction [V/m]
    - Ee     = Electric field in the east direction [V/m]
    - glat_e = Geographic latitude, accounting for ions [deg]
    - glon_e = Geographic longitude, accounting for ions [deg]
    '''

    shape_E = (Phi.shape[0]-2,Phi.shape[1]-2)
    En      = np.zeros(shape_E)
    Ee      = np.zeros(shape_E)
    lat_e   = np.zeros(shape_E)
    lon_e   = np.zeros(shape_E)
    wgs84   = nv.FrameE(name='WGS84') # for geographical position calculations; assume WGS84 ellipsoid

    for i in range(shape_E[0]):            # latitude index
        for j in range(shape_E[1]):        # longitude index
            lat_e[i,j] = glat[i+1,j+1]
            lon_e[i,j] = glon[i+1,j+1]
            # get 4 points on the surface
            point_A = wgs84.GeoPoint(latitude=glat[i,j+1]  ,longitude=glon[i,j+1]  ,z=z,degrees=True)
            point_B = wgs84.GeoPoint(latitude=glat[i+2,j+1],longitude=glon[i+2,j+1],z=z,degrees=True)
            point_C = wgs84.GeoPoint(latitude=glat[i+1,j]  ,longitude=glon[i+1,j]  ,z=z,degrees=True)
            point_D = wgs84.GeoPoint(latitude=glat[i+1,j+2],longitude=glon[i+1,j+2],z=z,degrees=True)
            # compute difference in distance between points
            # latitude
            delta_glat_m = point_A.delta_to(point_B).length
            # longitude
            delta_glon_m = point_C.delta_to(point_D).length
            # compute the E field as dPhi/dx
            En[i,j] = 0
            Ee[i,j] = 0
            if delta_glat_m!=0:
                En[i,j] = (Phi[i,j+1]-Phi[i+2,j+1])/delta_glat_m # direction of increasing AACGM lat
            if delta_glon_m!=0:
                Ee[i,j] = (Phi[i+1,j+2]-Phi[i+1,j])/delta_glon_m # direction of eastward AACGM lon

    return En,Ee,lat_e,lon_e
#_______________________________________________________________________________
def get_exb_drift_velocity(E:dict,B0:dict,Bu:dict,lat:np.ndarray,lon:np.ndarray,s=0.1): 
    '''Compute the horizontal E x B drift velocity component in the specified 
    direction at the given geographic latitude and longitude.  E is computed 
    from calculate_E_mix.  Assume B is vertical. Coordinate system is (east,north,up)  
    input: 
    - E     = Electric field, computed from calculate_E_mix [V/m] 
    - B0    = Magnitude of the magnetic field |B| [T] 
    - Bu    = The up component of the magnetic field [T]  
    - lat   = Geographic latitude [deg] 
    - lon   = Geographic longitude [deg] 
    - s     = Smoothing parameter for spline interpolation  
    output: 
    - v_ion_mag = Magnitude of the ion drift velocity [m/s] 
    - v_ion_dir = Direction of the ion drift velocity, defined as the angle between the east and north component [deg]
    - v_ion_e   = East component of the ion drift velocity [m/s]  
    - v_ion_n   = North component of the ion drift velocity [m/s]  
    '''

    # TODO: consider using the spline to calculate the derivatives directly 

    # convert to radians, flatten, and discard low-latitude boundary NaNs 
    theta   = np.pi - np.deg2rad(E['glat'][1:,:]).flatten() 
    phi     = np.deg2rad(E['glon'][1:,:]).flatten() 
    En_flat = E['En'][1:,:].flatten() 
    Ee_flat = E['Ee'][1:,:].flatten()

    theta_i         = np.pi - np.deg2rad(lat)
    phi_i           = np.deg2rad(lon)
    phi_i[phi_i<0] += 2.*np.pi 
    phi[phi<0]     += 2.*np.pi 

    # interpolation objects 
    En_interp_obj = scipy.interpolate.SmoothSphereBivariateSpline(theta,phi,En_flat,s=s)       
    Ee_interp_obj = scipy.interpolate.SmoothSphereBivariateSpline(theta,phi,Ee_flat,s=s)      

    th_min = np.min(theta) 
    th_max = np.max(theta) 
    ph_min = np.min(phi) 
    ph_max = np.max(phi) 

    # # get En, Ee at the interpolation point (theta_i,phi_i)
    # if( (theta_i.any()<th_min or theta_i.any()>th_max ) or 
    #     ( phi_i.any()<ph_min or phi_i.any()>ph_max ) ):
    #     # print("WARNING: glat range: th_min = {0:.5f}, th_max = {1:.5f}, th_i = {2:.5f}".format(th_min,th_max,theta_i[0])) 
    #     # print("WARNING: glon range: ph_min = {0:.5f}, ph_max = {1:.5f}, ph_i = {2:.5f}".format(ph_min,ph_max,phi_i[0]  ))
    #     print("WARNING: glat range: th_min = {0:.5f}, th_max = {1:.5f}, th_i out of range!".format(th_min,th_max)) 
    #     print("WARNING: glon range: ph_min = {0:.5f}, ph_max = {1:.5f}, ph_i out of range".format(ph_min,ph_max))
    #     print("Setting En = Ee = 0!") 
    #     En = 0
    #     Ee = 0 
    # else:
    #     En = En_interp_obj(theta_i,phi_i,grid=False)  
    #     Ee = Ee_interp_obj(theta_i,phi_i,grid=False) 
    En = En_interp_obj(theta_i,phi_i,grid=False)  
    Ee = Ee_interp_obj(theta_i,phi_i,grid=False) 

    # calculate the ion drift vector
    v_ion_e =       En*Bu/(B0*B0) 
    v_ion_n = (-1.)*Ee*Bu/(B0*B0)

    # get magnitude and direction 
    v_ion_mag = np.squeeze( np.sqrt( np.power(v_ion_e,2.) + np.power(v_ion_n,2.) ) )
    if(v_ion_n.any()!=0):
        # degrees East of North 
        v_ion_dir = np.squeeze( np.rad2deg( np.arctan(v_ion_e/v_ion_n) ) ) 
    else:
        v_ion_dir = 0  

    return v_ion_mag,v_ion_dir,v_ion_e,v_ion_n  
