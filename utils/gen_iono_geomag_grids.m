function [iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    Bx, By, Bz, geomag_grid_parms] = gen_iono_geomag_grids(...
    hts, lats, lons, time, R12, B_ht_inc, B_lat_inc, B_lon_inc, iri_options)

arguments 
    hts double
    lats double
    lons double
    time double
    R12 double
    B_ht_inc double
    B_lat_inc double
    B_lon_inc double
    iri_options struct = []
end


%% generate ionospheric, geomagnetic and irregularity grids
% [iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
%     Bx, By, Bz, geomag_grid_parms] = ...
%     gen_iono_geomag_grids(hts, lats, lons, time, R12, B_ht_inc, B_lat_inc, B_lon_inc)


%% Specify time in PHaRLaP format
year = str2double(datestr(time, 'yyyy'));
month = str2double(datestr(time, 'mm'));
day  = str2double(datestr(time, 'dd'));
hour = str2double(datestr(time, 'HH'));
minute = str2double(datestr(time, 'MM'));
UT = [year month day hour minute];

% geomag
doppler_flag = 0;

ht_start = hts(1);          % start height for ionospheric grid (km)
ht_inc = unique(diff(hts));             % height increment (km)
num_ht = length(hts);
lat_start = lats(1);
lat_inc = uniquetol(diff(lats), 0.001);
num_lat = length(lats);
lon_start = lons(1);
lon_inc = uniquetol(diff(lons), 0.001);
lon_inc = abs(lon_inc(1));
num_lon = length(lons);


assert(num_ht < 402)
assert(num_lat < 102)
assert(num_lon < 102)

B_ht_start = ht_start;          % start height for geomagnetic grid (km)
B_num_ht = ceil(num_ht .* ht_inc ./ B_ht_inc);
B_lat_start = lat_start;
B_num_lat = ceil(num_lat .* lat_inc ./ B_lat_inc);
B_lon_start = lon_start;
B_num_lon = ceil(num_lon .* lon_inc ./ B_lon_inc);
if isempty(lon_inc)
    B_num_lon = 1;
end
assert(B_num_ht < 202)
assert(B_num_lat < 102)
assert(B_num_lon < 102)

iono_grid_parms = [lat_start, lat_inc, num_lat, lon_start, lon_inc, num_lon, ...
    ht_start, ht_inc, num_ht, ];

geomag_grid_parms = [B_lat_start, B_lat_inc, B_num_lat, B_lon_start, ...
    B_lon_inc, B_num_lon, B_ht_start, B_ht_inc, B_num_ht];

tic
fprintf('Generating ionospheric and geomag grids... ')

if isempty(iri_options)
    [iono_pf_grid, iono_pf_grid_5, collision_freq, Bx, By, Bz] = ...
        gen_iono_grid_3d(UT, R12, iono_grid_parms, geomag_grid_parms, ...
        doppler_flag);

else
    [iono_pf_grid, iono_pf_grid_5, collision_freq, Bx, By, Bz] = ...
        gen_iono_grid_3d(UT, R12, iono_grid_parms, geomag_grid_parms, ...
        doppler_flag,  'iri2016', iri_options);

end

% [iono_pf_grid, iono_pf_grid_5, collision_freq, Bx, By, Bz] = ...
%     gen_iono_grid_3d(UT, R12, iono_grid_parms, geomag_grid_parms, ...
%     doppler_flag);
toc
fprintf('\n')

% convert plasma frequency grid to  electron density in electrons/cm^3
iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;


