from nc_utils import ncread_vars
import matplotlib.pyplot as plt
import os
data = ncread_vars(os.path.expanduser('~/Downloads/test_2.nc'))
plt.contourf(data['lon'], data['lat'], data['hipc'][50, :, :])
plt.show()
