%% inputs

time = datenum(2020, 7, 14);
freq = 162; % 2:0.2:20; % 20;
foEs = 20;

alts = 80:2:400;
sat_alt = 890;

txlat = 36;
txlon = -58;
rxlat = 36;
rxlon = -76;
alt = 0;

% grid
lats = 34:38;
lons = -80:2:46;

elvarr = [2:30];
azarr = -100:-80;
R12 = 100;

OX_mode = 0;

kp = 3;
maxdist = 1E3;  % meters from homing
tol = [1e-7 0.01 25];       % ODE solver tolerance and min max stepsizes


%% Specify time in PHaRLaP format
year = str2double(datestr(time, 'yyyy'));
month = str2double(datestr(time, 'mm'));
day  = str2double(datestr(time, 'dd'));
hour = str2double(datestr(time, 'HH'));
minute = str2double(datestr(time, 'MM'));
UT = [year month day hour minute];


%% Load ionosphere
iri_options = [];
iri_options.foE = foEs;
iri_options.profile_type = 'iri2016';
[iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    Bx, By, Bz, geomag_grid_parms] = ...
    gen_iono_geomag_grids(alts, lats, lons, UT, R12, iri_options);
iono_pf_grid = sqrt(80.6 * iono_en_grid ./ 1E6);


%% Raytrace
txloc = [txlat, txlon, alt];
rxloc = [rxlat, rxlon, alt];
% homed_ray = ...
%     raytrace_itsi(freq, OX_mode, txloc, rxloc, iono_en_grid, iono_en_grid_5, ...
%     collision_freq, iono_grid_parms, Bx, By, Bz, geomag_grid_parms, ...
%     elvarr, azarr, maxdist, tol, 1, 0);


origin_lat = txloc(1);
origin_long = txloc(2);
origin_ht = txloc(3);

nhops = 1; % number of hops

rays = gs_raytrace(elvarr, azarr, freq, nhops, OX_mode, ...
    origin_lat, origin_long, origin_ht, iono_en_grid, iono_en_grid_5, ...
    collision_freq, iono_grid_parms, Bx, By, Bz, geomag_grid_parms, tol, refractive_ind);

plot_rays(txloc, rxloc, rays)



























