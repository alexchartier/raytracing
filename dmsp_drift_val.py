import numpy as np
import h5py
import os
import datetime as dt

"""
Validate potential maps using DMSP cross-track drift data
"""

def main():
    in_fname = os.path.expanduser('~/Downloads/dms_20190301_15s1.001.hdf5')
    preproc_dmsp(in_fname)

def preproc_dmsp(in_fname):
    hf = h5py.File(in_fname, 'r')
    vals = {}
    headerdata = hf['Metadata']['Data Parameters'][...]
    header = [h[0].decode('UTF-8') for h in headerdata]
    assert 'NE' in header, 'NE not in header'
    data = hf['Data']['Table Layout'][...]  # 1Hz DMSP data
   
    breakpoint()
    for h in header:
        vals[h] = np.array([])
    vals['time'] = np.array([])

    for dind, d in enumerate(data):
        for hind, h in enumerate(header):
            vals[h] = np.append(vals[h], d[hind])
        time = dt.datetime(int(vals['YEAR'][dind]), int(vals['MONTH'][dind]), int(vals['DAY'][dind]), \
                           int(vals['HOUR'][dind]), int(vals['MIN'][dind]), int(vals['SEC'][dind]))

        vals['time'] = np.append(vals['time'], time)


if __name__ == '__main__':
    main()
