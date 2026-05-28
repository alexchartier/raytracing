function [iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    Bx, By, Bz, geomag_grid_parms] = gen_grid_parms(sami_t)
%% generate ionospheric, geomagnetic and irregularity grids from SAMI3

% unpack
hts = sami_t.alt;
lats = sami_t.lat;
lons = sami_t.lon;
iono_en_grid = sami_t.dene;
iono_en_grid_5 = sami_t.dene;


ht_start = hts(1);          % start height for ionospheric grid (km)
ht_inc = unique(diff(hts));             % height increment (km)
num_ht = length(hts);
lat_start = lats(1);
lat_inc = uniquetol(diff(lats), 0.001);
num_lat = length(lats);
lon_start = lons(1);
lon_inc = uniquetol(diff(lons), 0.001);
num_lon = length(lons);

% geomag
B_ht_inc = ht_inc;                  % height increment (km)
B_lat_inc = lat_inc;
B_lon_inc = lon_inc;

B_ht_start = ht_start;          % start height for geomagnetic grid (km)
B_num_ht = ceil(num_ht .* ht_inc ./ B_ht_inc);
B_lat_start = lat_start;
B_num_lat = ceil(num_lat .* lat_inc ./ B_lat_inc);
B_lon_start = lon_start;
B_num_lon = ceil(num_lon .* lon_inc ./ B_lon_inc);

iono_grid_parms = [lat_start, lat_inc, num_lat, lon_start, lon_inc, num_lon, ...
    ht_start, ht_inc, num_ht, ];

geomag_grid_parms = [B_lat_start, B_lat_inc, B_num_lat, B_lon_start, ...
    B_lon_inc, B_num_lon, B_ht_start, B_ht_inc, B_num_ht];

Bx = zeros(size(iono_en_grid));
By = zeros(size(iono_en_grid));
Bz = zeros(size(iono_en_grid));
collision_freq = zeros(size(iono_en_grid));


