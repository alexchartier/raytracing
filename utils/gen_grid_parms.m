function [iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    Bx, By, Bz, geomag_grid_parms] = gen_grid_parms(model)
%% generate ionospheric, geomagnetic and irregularity grids from SAMI3
% model.dene = alts/lats/lons order
% model.lat - degrees
% model.lon - degrees
% model.alt - km
% all regularly spaced


% %% Check for 0/360 issues
% 
%     loni = model.lon < 360;
%     model.dene = model.dene(:, :, loni);
%     model.lon = model.lon(:, :, loni);
% end

%% unpack
hts = model.alt;
lats = model.lat;


if min(model.lon) <= 0 && max(model.lon) >= 360
    loni = model.lon < 360;
    lons = model.lon(loni);
    dene = model.dene(:, :, loni);
    li = lons >= 180;
    lons = [lons(li) - 360, lons(~li)];

    iono_en_grid = permute(dene, [2, 3, 1]);
    iono_en_grid = [iono_en_grid(:, li, :), iono_en_grid(:, ~li, :)];
    iono_en_grid = [iono_en_grid, iono_en_grid(:, 1, :)];
    lons = [lons, lons(1) + 360];
else

    lons = model.lon;
    li = model.lon >= 180;
    lons = [lons(li) - 360, lons(~li)];

    iono_en_grid = permute(model.dene, [2, 3, 1]);
    iono_en_grid = [iono_en_grid(:, li, :), iono_en_grid(:, ~li, :)];
end

iono_en_grid_5 = iono_en_grid;

ht_start = hts(1);          % start height for ionospheric grid (km)
ht_diff = diff(hts);
if isempty(ht_diff)
    ht_inc = 0;
else
    ht_inc = ht_diff(1);
end
num_ht = length(hts);
lat_start = lats(1);
lat_diff = diff(lats);
if isempty(lat_diff)
    lat_inc = 0;
else
    lat_inc = lat_diff(1);
end
num_lat = length(lats);
lon_start = lons(1);
lon_diff = diff(lons);
if isempty(lon_diff)
    lon_inc = 0;
else
    lon_inc = lon_diff(1);
end
num_lon = length(lons);

iono_grid_parms = [lat_start, lat_inc, num_lat, lon_start, lon_inc, num_lon, ...
    ht_start, ht_inc, num_ht];

%% geomag
if length(hts) <= 101

B_ht_inc = ht_inc;                  % height increment (km)
else    
    B_ht_inc = (max(hts) - min(hts)) / 50;
    hts = min(hts):B_ht_inc:max(hts);
    
end
B_lat_inc = lat_inc;
B_lon_inc = lon_inc;

B_ht_start = ht_start;          % start height for geomagnetic grid (km)
B_num_ht = length(hts);
B_lat_start = lat_start;
B_num_lat = ceil(num_lat .* lat_inc ./ B_lat_inc);
B_lon_start = lon_start;
B_num_lon = ceil(num_lon .* lon_inc ./ B_lon_inc);

geomag_grid_parms = [B_lat_start, B_lat_inc, B_num_lat, B_lon_start, ...
    B_lon_inc, B_num_lon, B_ht_start, B_ht_inc, B_num_ht];


[Bx, By, Bz] = gen_bfield(lats, lons, hts, year(model.time));

%% Plot the B-field



%% Collisions
collision_freq = zeros(size(iono_en_grid));
































