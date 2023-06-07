%%

fname = '~/Downloads/test_2.nc';

ncdisp(fname)
hipc = ncread(fname, 'hipc');
hihc = ncread(fname, 'hihc');
lats =  ncread(fname, 'lat');
lons =  ncread(fname, 'lon');

cmax = 30;
minlat = 50;
subplot(2, 1, 1)
contourf(lons, lats(lats > minlat), hipc(:, (lats > minlat), 100)', 60)

clim([0, cmax])
title('hipc')

subplot(2, 1, 2)
contourf(lons, lats(lats > minlat), hihc(:, (lats > minlat), 100)', 60)
clim([0, cmax])
title('hihc')
