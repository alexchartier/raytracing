import matplotlib.pyplot as plt
import nc_utils
import datetime as dt
import numpy as np
from scipy.interpolate import griddata
from dmsp_drift_val import get_rad_theta


def main(
    time = dt.datetime(2019, 1, 2, 10, 20),
    in_fname = '~/data/sami3/2019/sami3_regulargrid_conductances_2019Mar02.nc',
):

    ti = 0
    data = nc_utils.ncread_vars(in_fname)

    pts = (data['cond_lat'].flatten(), data['cond_lon'].flatten())
    vals = data['hipc'][ti, :, :].flatten()
    lats_i = np.arange(0, 91)
    lons_i = np.arange(0, 361)
    xi, yi = np.meshgrid(lats_i, lons_i)
    vals_i = griddata(pts, vals, (xi.flatten(), yi.flatten()), )
    vals_ig = vals_i.reshape((len(lons_i), len(lats_i)))



    fig, ax = plt.subplots(1, 1,subplot_kw={'projection': 'polar'})

    latind = lats_i > 40
    rad, theta = get_rad_theta(lats_i[latind], lons_i)
    plt.contourf(theta, rad, vals_ig[:, latind].T)
    plt.show()


if __name__ == '__main__':
    main()
