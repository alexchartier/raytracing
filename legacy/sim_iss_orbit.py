import ephem
import numpy as np
import datetime as dt
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

name = "ISS (ZARYA)"
line1 = "1 25544U 98067A   21322.45522487  .00001781  00000+0  41042-4 0  9998"
line2 = "2 25544  51.6445 299.9201 0004682 212.2898 289.6236 15.48591869312504"

tle_rec = ephem.readtle(name, line1, line2)

stime = dt.datetime(2021, 11, 18, 12, 0)
etime = dt.datetime(2021, 11, 18, 13, 40)
tstep = dt.timedelta(minutes=1)

time = stime
ax = plt.axes(projection=ccrs.PlateCarree())
while time < etime:
    tle_rec.compute(time)
    print(time, np.rad2deg(tle_rec.sublat), np.rad2deg(tle_rec.sublong))
    time += tstep
    plt.plot(np.rad2deg(tle_rec.sublong), np.rad2deg(tle_rec.sublat), '.', color='red')

ax.coastlines()

plt.show()


