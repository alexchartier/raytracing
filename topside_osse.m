%% topside osse
% Simulate 2-satellite sounding


%% input params

UT = [2000 9 21 0 0];           % UT - year, month, day, hour, minute
speed_of_light = 2.99792458e8;
R12 = 100;
%freqs = [2:8, 8.1:0.1:12, 13:14];   % frequency (MHz)
freqs = [10:0.2:13];   % frequency (MHz)
rxloc = [0.0, 130.0, 600.0];  % transmitter location
txloc = [-10.0, 130, 600];  % receiver location
tol = [1e-8 0.01 25];       % ODE solver tolerance and min max stepsizes

doppler_flag = 0;               % interested in Doppler shift

fprintf( ['\n' ...
    'Topside example of 3D magneto-ionic numerical raytracing for a WGS84 ellipsoidal' ...
    ' Earth\n\n'])
% 
% % Define iono and geomag grids  * make sure these go into 10km and 1 degree
% hts = 200:1:900;
% lats = -12:0.2:-3;
% lons = 128:0.2:132;


%% load/generate ionospheric, geomagnetic and irregularity grids
% [iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
%     Bx, By, Bz, geomag_grid_parms] = gen_iono_geomag_grids(hts, lats, lons, UT, R12);


% generate ionospheric, geomagnetic and irregularity grids
%
ht_start = 200;          % start height for ionospheric grid (km)
ht_inc = 1;             % height increment (km)
num_ht = 401;         
hts = ht_start:ht_inc:(ht_start + ht_inc * (num_ht - 1));
lat_start = -20.0;
lat_inc = 0.3;
num_lat = 101.0;
lats = lat_start:lat_inc:(lat_start + lat_inc * (num_lat - 1));
lon_start= 128.0;
lon_inc = 1.0;
num_lon = 5.0;
iono_grid_parms = [lat_start, lat_inc, num_lat, lon_start, lon_inc, num_lon, ...
      ht_start, ht_inc, num_ht, ];

B_ht_start = ht_start;          % start height for geomagnetic grid (km)
B_ht_inc = 5;                  % height increment (km)
B_num_ht = ceil(num_ht .* ht_inc ./ B_ht_inc);
b_hts = B_ht_start:B_ht_inc:(B_ht_start + B_ht_inc * (B_num_ht - 1));
B_lat_start = lat_start;
B_lat_inc = 1.0;
B_num_lat = ceil(num_lat .* lat_inc ./ B_lat_inc);
B_lon_start = lon_start;
B_lon_inc = 1.0;
B_num_lon = ceil(num_lon .* lon_inc ./ B_lon_inc); 
geomag_grid_parms = [B_lat_start, B_lat_inc, B_num_lat, B_lon_start, ...
      B_lon_inc, B_num_lon, B_ht_start, B_ht_inc, B_num_ht];

tic
fprintf('Generating ionospheric and geomag grids... ')
[iono_pf_grid, iono_pf_grid_5, collision_freq, Bx, By, Bz] = ...
    gen_iono_grid_3d(UT, R12, iono_grid_parms, ...
                     geomag_grid_parms, doppler_flag);
toc
fprintf('\n')

% convert plasma frequency grid to  electron density in electrons/cm^3
iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;



%% homing raytrace over all freqs
rays = [];

az = azimuth(txloc(1), txloc(2), rxloc(1), rxloc(2));
elvarr = -45:0.5:0;
azarr = az; % - 2:2:az + 2;
maxdist = 50E3;
freqs = 45;
for f = 1%1:length(freqs)
    fprintf('Freq: %1.1f MHz\n', freqs(f)) 
    rays_f = raytrace_bistatic(freqs(f), txloc, rxloc, iono_en_grid, iono_en_grid_5, ...
    collision_freq, iono_grid_parms, Bx, By, Bz, geomag_grid_parms, elvarr, azarr, maxdist, tol);
    rays = [rays, rays_f];
end


%%
dene = squeeze(iono_en_grid(:, 5, :))';
hold on
contourf(lats, hts, dene)
for r = 1:length(rays)
    plot(rays(r).lat, rays(r).height, 'w')
end

plot(txloc(1), txloc(3), 'r.', 'markersize', 20)
plot(rxloc(1), rxloc(3), 'g.', 'markersize', 20)
% plot_rays(txloc, rxloc, rays)





