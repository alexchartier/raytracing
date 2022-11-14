"""
Convert a sami3 output file to a netCDF electron density file to use for raytracing in HFSIM
"""
import sys
from datetime import datetime
import numpy as np
import os
import errno
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import nc_utils
from scipy.interpolate import SmoothSphereBivariateSpline

LAT_DIMENSION = 90
LON_DIMENSION = 180
ALT_DIMENSION = 150

def main(in_fname):
    sami = sami_to_icd(in_fname)

    # define the write variables
    dim_map = get_dim_map()
    dim_defs = def_dims(sami, dim_map)
    out_vars = get_vars(sami, dim_map)

    nc_utils.write_nc(out_fname, out_vars, dim_defs)


def sami_to_icd(sami_fname):
    sami = nc_utils.ncread_vars(sami_fname)
    for k, v in sami.items():
        sami[k] = v.astype('float64')

    for ia, alt in enumerate(sami['alt0']):
        lats = np.deg2rad(sami['lat0'][ia, :, :]).ravel() + np.pi / 2
        lons = np.deg2rad(sami['lon0'][ia, :, :]).ravel()
        data = sami['dene0'][0, ia, :, :].T.ravel()
        breakpoint()
        lut = SmoothSphereBivariateSpline(lats, lons, data, s=3.5)
        for it, time in enumerate(sami['time']):
            data = sami['dene0'][it, ia, :, :]
            lut.F = data 

            # TODO add pi/2 to interp lats

    # remove singleton dimension
    sami['dene'] = np.squeeze(sami['dene'])
    
    # Get the alt, lat, lon array for the specified time
    hourDecimal = time.hour + time.minute/60
    hours = np.asarray(sami['hrut'])
    timeIdx = (np.abs(hours - hourDecimal)).argmin()
    sami['dene'] = sami['dene'][timeIdx, :, :, :]

    # remove duplicated 180 in longitude
    sami['lon'] = sami['lon'][1:]
    sami['dene'] = sami['dene'][:, :, 1:]

    # Straighten out the negatives (-180 to 180 vs 0 to 360)
    neglon_idx = sami['lon'] < 0
    sami['dene'] = np.concatenate((sami['dene'][:, :, ~neglon_idx], sami['dene'][:, :, neglon_idx]), axis=2)
    sami['lon'][neglon_idx] += 360
    sami['lon'] = np.concatenate((sami['lon'][~neglon_idx], sami['lon'][neglon_idx]), axis=0)

    sami['dene'] *= 1E6  # el. cm3 -> el. m3

    # sami['dene'][sami['dene'] == 0] = 33333E6  # remove zeros
   
    # shift from alt/lat/lon to lat/lon/alt 
    sami['dene'] = np.squeeze(sami['dene']).transpose((1, 2, 0))

    # fix dimensions to 50x50x100  (to match with the HFSIM required dimensions)
    sami['lon'] = sami['lon'][-50:]
    sami['lat'] = sami['lat'][-50:]
    if len(sami['alt']) < ALT_DIMENSION:
        # Interpolate to extend to the required number of altitude data points
        startingArray = sami['alt']
        newArrayLength = ALT_DIMENSION
        delta = (len(startingArray)-1) / (newArrayLength-1)
        sami['alt'] = [interpolate(startingArray, i*delta) for i in range(newArrayLength)]

        newElectronDensityArray = np.zeros((LAT_DIMENSION, LON_DIMENSION, ALT_DIMENSION))

        # Interpolate the electron density data for the new altitudes
        for i in range(-LAT_DIMENSION,-1):
            for j in range(-LON_DIMENSION,-1):
                # Get the current electron desity profile for the given lat/lon
                startingArray = sami['dene'][i, j, :]
                newArrayLength = ALT_DIMENSION
                delta = (len(startingArray)-1) / (newArrayLength-1)
                # Interpolate the data to the desired length
                newElectronDensities = [interpolate(startingArray, i*delta) for i in range(newArrayLength)]
                newElectronDensityArray[i,j] = newElectronDensities
        
        sami['dene'] = newElectronDensityArray

    else:
        sami['alt'] = sami['alt'][:100]
        sami['dene'] = sami['dene'][-50:, -50:, :100]

    # Loop through each altitude and smooth over any zero values
    for i in range(len(sami['alt'])):
        altGrid = sami['dene'][:, :, i]
        min = np.min(altGrid[np.nonzero(altGrid)])
        max = np.max(altGrid[np.nonzero(altGrid)])
        for j in range(len(sami['lat'])):
            altLatRow = np.array(altGrid[j,:])

            # If the row has no zeros, skip the row
            # If the row is all zeros, skip it because you can't interpolate
            #   all zeros
            if 0 not in altLatRow or np.all((altLatRow == 0)):
                continue

            k = np.arange(len(altLatRow))
            idx = np.nonzero(altLatRow)
            f = interp1d(k[idx],altLatRow[idx],kind='linear',fill_value="extrapolate")
            #f = interp1d(k[idx],altLatRow[idx],fill_value=(min), bounds_error=False)
            altLatRowNew = f(k)

            altGrid[j,:] = altLatRowNew

        # If there are still zeros because a row was skipped, redo interpolation
        # by column 
        if 0 in altGrid:
            for j in range(len(sami['lat'])):
                altLonRow = np.array(altGrid[:,j])

                # If the row has no zeros, skip the row
                if 0 not in altLonRow:
                    continue

                k = np.arange(len(altLonRow))
                idx = np.nonzero(altLonRow)
                f = interp1d(k[idx],altLonRow[idx],kind='nearest',fill_value="extrapolate")
                altLonRowNew = f(k)

                altGrid[:,j] = altLonRowNew 

        sami['dene'][:, :, i] = altGrid

            # # f = interp2d(sami['lat'], sami['lon'], altGrid, kind='cubic')
            # f = interp2d(np.nonzero(altGrid)[0], np.nonzero(altGrid)[1], altGrid[np.nonzero(altGrid)])
            # interpolatedAltGrid = f(sami['lat'], sami['lon'])
            # sami['dene'][:,:,i] = interpolatedAltGrid


    return sami



if __name__ == '__main__':

    """
    args = sys.argv

    assert len(args) == 4, 'Should have 3 arguments, e.g.:\n' + \
        'python3 sami3_to_hfsim.py 2022,03,12,16,15 /path/to/input/file.nc ' + \
        '/path/to/output/file.nc'
    
    time = datetime.strptime(args[1], '%Y,%m,%d,%H,%M')
    inputDir = args[2]
    outputFilename = args[3]

    """

    in_fname = 'sami3_NEXT.nc'
    main(in_fname)
