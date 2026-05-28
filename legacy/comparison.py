import numpy as np
import matplotlib.pyplot as plt
import nvector as nv


def main():
    lat = np.arange(60, 75)
    lon = np.zeros(lat.shape)
    alt = np.ones(lat.shape) * 800E3
    z = -alt
    deg = True

    # convert geodetic (lat,lon,z) to ECEF (x,y,z) using nvector
    wgs84 = nv.FrameE(name='WGS84')
    N = len(lat) - 1
    XYZ = np.zeros((N, 3))
    UVW = np.zeros((N, 3))
    for i in range(N):
        p = wgs84.GeoPoint(
            latitude=lat[i], longitude=lon[i], z=z[i], degrees=deg)
        p2 = wgs84.GeoPoint(
            latitude=lat[i+1], longitude=lon[i+1], z=z[i+1], degrees=deg)

        sat_brng = p.delta_to(p2) 
        p_ecef = p.to_ecef_vector()
        sat_brng_ecef = sat_brng.to_ecef_vector()

        #TODO calc perpendicular bearing
        n_EB_E = sat_brng
        frame_B = nv.FrameB(sat_brng, yaw=90, degrees=True)  # body frame, looking perpendicular right to satellite
        p_BC_B = frame_B.Pvector(np.r_[1000, 0, 0].reshape(-1, 1))
        p_BC_E = p_BC_B.to_ecef_vector()
        p_EB_E = 

        perp_brng = p.delta_to(perp_pt)

        XYZ[i, :] = p_ecef.pvector[:].flatten()
        UVW[i, :] = sat_brng_ecef.pvector[:].flatten()

        breakpoint()

    # plot the satellite bearings
    ax = plt.figure().add_subplot(projection='3d')
    ax.quiver(*XYZ.T, *UVW.T)
    plt.show()


if __name__ == '__main__':
    main()
