function [iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    Bx, By, Bz, geomag_grid_parms] = gen_iono_geomag_grids(D)

%% generate ionospheric, geomagnetic and irregularity grids
% % Define iono and geomag grids  * make sure these go into 10km and 1 degree

%% iono

doppler_flag = 0;

ht_start = D.alt(1);          % start height for ionospheric grid (km)
ht_inc = unique(diff(D.alt));             % height increment (km)
num_ht = length(D.alt);
lat_start = D.lat(1);
lat_inc = uniquetol(diff(D.lat), 0.001);
num_lat = length(D.lat);
lon_start= D.lon(1);
lon_inc = uniquetol(diff(D.lon), 0.001);
num_lon = length(D.lon);

assert(num_ht < 402, 'raytrace_3d cannot support more than 401 heights')
assert(num_lat < 702, 'raytrace_3d cannot support more than 701 lats')
assert(num_lon < 702, 'raytrace_3d cannot support more than 701 lons')

iono_grid_parms = [lat_start, lat_inc, num_lat, lon_start, lon_inc, num_lon, ...
    ht_start, ht_inc, num_ht, ];

geomag_grid_parms = iono_grid_parms;

iono_en_grid = permute(D.dene, [3, 2, 1]);
iono_en_grid_5 = iono_en_grid;
collision_freq = zeros(size(iono_en_grid));
Bx = zeros(size(iono_en_grid));
By = zeros(size(iono_en_grid));
Bz = zeros(size(iono_en_grid));
iono_pf_grid = sqrt(80.6 * iono_en_grid ./ 1E6);

assert(size(Bx, 1) <= 101, 'Bx, By, Bz max 101 heights')
assert(size(Bx, 2) <= 101, 'Bx, By, Bz max 101 lats')
assert(size(Bx, 3) <= 201,  'Bx, By, Bz max 201 lons')





























