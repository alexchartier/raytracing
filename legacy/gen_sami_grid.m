% function [iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
%     Bx, By, Bz, geomag_grid_parms] = ...
%     gen_sami_grid(hts, lats, lons, sami_fname, UT, R12, varargin);

%% generate ionospheric, geomagnetic and irregularity grids
% % Define iono and geomag grids  * make sure these go into 10km and 1 degree
hts = 100:2:900;
lats = -12:0.5:-3;
lons = 128:0.5:132;
sami_fn = 'sami3_NEXT.nc'; 
UT = [2000 9 21 0 0];           % UT - year, month, day, hour, minute
R12 = 100;
% [iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, Bx, By, Bz, geomag_grid_parms] = ...
%     gen_sami_grid(hts, lats, lons, sami_fname, UT, R12);

time_hr = UT(3);

%% load SAMI3

D = load_sami(sami_fn);

%% Interpolate to regular lat/lon
[min_t, ti] = min(abs(D.time - time_hr));
assert(min_t < 0.2, 'Times too far away')

dene1 = ones(length(D.alt0), length(lats), length(lons));

for a = 1:length(D.alt0)
    F = scatteredInterpolant(D.lat0(a, :, :), D.lon0(a, :, :), D.dene0(a, :, :, ti));
    dene1(a, :, :) = F(lats, lons);

end

%Tind = 


%% iono grid 

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