import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation as R
from spacepy import pycdf
import os
import xarray
import nc_utils

"""
=====================================================
Rotating E-field from XYZ to North-East-Centre (NEC):
Rotation matrix based on the sat. velocity in NEC

=====================================================
"""


def main(
    swarm_fn = '~/raytracing/data/swarm/SW_EXPT_EFIA_TCT02_20190302T000000_20190302T121504_0302.cdf',
    nc_fname_fmt = 'data/swarm_proc/%Y%b%d_swarm.nc',
):
    df = read_swarm(swarm_fn)
    nc_fname = df.index[0].strftime(nc_fname_fmt)
    nc_utils.write_netcdf_from_df(df, nc_fname)   


def read_swarm(swarm_fn, qualflag=4):
    """ reads swarm file into a dataframe, with NEC ion drift conversion 

         Bit 0 (least significant): 𝑣𝑖,𝑥,𝐻 (along-track component from the horizontal sensor)
         Bit 1: 𝑣𝑖,𝑥,𝑉 (along-track component from the vertical sensor)
         Bit 2: 𝑣𝑖,𝑦 (to the right, observer facing forward)
         Bit 3: 𝑣𝑖,𝑧 (downward)
    """
    vi_nec = []
    TII_data = pd.DataFrame.from_dict(nc_utils.load_cdf(swarm_fn)).set_index('Timestamp')
    TII_data = TII_data[TII_data['Quality_flags'] == qualflag]    
    assert max(TII_data['Quality_flags']) == 4, 'Need to check assumptions re. qual flags'

    for i in range(len(TII_data)):
        mat = R.align_vectors(np.array([[TII_data['VsatN'][i], TII_data['VsatE'][i], TII_data['VsatC'][i]]]), np.array([[1, 0, 0]]))
        mat = mat[0].as_matrix()     
        #E_nec.append(mat.dot([TII_data['Ex'][i], TII_data['Ey'][i], TII_data['Ez'][i]]))
        vi_nec.append(mat.dot([0, TII_data['Viy'][i], 0]))

    vi_nec = pd.DataFrame(vi_nec, columns=['Vn','Ve','Vc']).set_index(TII_data.index)
    
    for k in 'Latitude', 'Longitude', 'Radius', 'Viy', 'Viy_error':
        vi_nec[k] = TII_data[k]
    
    return vi_nec



if __name__ == '__main__':
    main()






