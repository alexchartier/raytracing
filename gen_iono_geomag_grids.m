function [iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    Bx, By, Bz, geomag_grid_parms] = gen_iono_geomag_grids(hts, lats, lons, UT, R12)
%% % generate ionospheric, geomagnetic and irregularity grids
% % Define iono and geomag grids  * make sure these go into 10km and 1 degree
% hts = 100:2:900;
% lats = -12:0.5:-3;
% lons = 128:0.5:132;
% UT = [2000 9 21 0 0];           % UT - year, month, day, hour, minute
% R12 = 100;
%[iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, Bx, By, Bz, geomag_grid_parms] = gen_iono_geomag_grids(hts, lats, lons, UT, R12)


%% geomag

doppler_flag = 0;

ht_start = hts(1);          % start height for ionospheric grid (km)
ht_inc = unique(diff(hts));             % height increment (km)
num_ht = length(hts);
lat_start = lats(1);
lat_inc = uniquetol(diff(lats), 0.001);
num_lat = length(lats);
lon_start= lons(1);
lon_inc = uniquetol(diff(lons), 0.001);
num_lon = length(lons);

assert(num_ht < 402, 'raytrace_3d cannot support more than 401 heights')
assert(num_lat < 702, 'raytrace_3d cannot support more than 701 lats')
assert(num_lon < 702, 'raytrace_3d cannot support more than 701 lons')

iono_grid_parms = [lat_start, lat_inc, num_lat, lon_start, lon_inc, num_lon, ...
    ht_start, ht_inc, num_ht, ];


B_ht_start = ht_start;          % start height for geomagnetic grid (km)
B_ht_inc = 10;                  % height increment (km)
B_num_ht = ceil(num_ht .* ht_inc ./ B_ht_inc);
B_lat_start = lat_start;

if rem(length(lats) - 1, 2) == 0
    B_lat_inc = 2.0;
elseif rem(length(lats) - 1, 3) == 0
    B_lat_inc = 3.0;
elseif rem(length(lats) -1, 5) == 0
    B_lat_inc = 5.0;
else
    fprintf('Pick a number of lats that divides into 2, 3,  or 5')
end

B_num_lat = ceil(num_lat .* lat_inc ./ B_lat_inc);
B_lon_start = lon_start;

if rem(length(lons) - 1, 2) == 0
    B_lon_inc = 2.0;
elseif rem(length(lons) -1, 3) == 0
    B_lon_inc = 3.0;
elseif rem(length(lons) - 1, 5) == 0
    B_lon_inc = 5.0;
else
    fprintf('Pick a number of lons that divides into 2, 3,  or 5')
end
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



assert(size(Bx, 1) <= 101, 'Bx, By, Bz max 101 heights')

assert(size(Bx, 2) <= 101, 'Bx, By, Bz max 101 lats')
assert(size(Bx, 3) <= 201,  'Bx, By, Bz max 201 lons')





























