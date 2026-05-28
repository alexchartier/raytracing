%%

fname = '~/Downloads/test_2.nc';
fname_in = '/Users/chartat1/pymix/data/sami3_conductance/sami3_NEXT.nc';
ncdisp(fname_in)


%%
hipc = ncread(fname, 'hipc');
hihc = ncread(fname, 'hihc');
lats =  ncread(fname, 'lat');
lons =  ncread(fname, 'lon');
% 
lats_int = -89:1:89;
lons_int = 0:2:360;
lat_in =  ncread(fname_in, 'lat');
lon_in =  ncread(fname_in, 'lon');
hipcS_in = ncread(fname_in, 'hipcS');
idx = [1:5:30, 31:80];
lon_in = double(squeeze(lon_in(1, idx, :)));
lat_in = double(squeeze(lat_in(1, idx, :)));
hipcS_in = double(squeeze(hipcS_in(idx, :, 100)));


% %% Try out some new interpolation
[lon2d, lat2d] = meshgrid(lons_int, lats_int);
F = scatteredInterpolant(lon_in(:), lat_in(:), hipcS_in(:));
hipc_int = F(lon2d, lat2d);


%
clf
cmax = 550;
minlat = 50;
subplot(2, 1, 1)
hold on
n = 1;
%h  = pcolor(lons(1:n:end), lats(1:n:end), hipc(1:n:end, 1:n:end, 100)');
[~, h] = contourf(lons, lats, squeeze(hipc(:, :, 100))');
% [~, h]  = contourf(lons_int(1:end), lats_int, squeeze(hipc_int(:, 1:end)));
set(h, 'EdgeColor', 'none')
%plot(lon_in(:), lat_in(:), '.k')
hold off

%clim([0, cmax])
title('hipc')

subplot(2, 1, 2)
h = pcolor(lons_int, lats_int, hipc_int);
set(h, 'EdgeColor', 'none')
% clim([0, cmax])