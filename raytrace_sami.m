%% inputs
time = datetime(2014, 5, 23, 6, 0, 0);
sami_fn = 'test.nc';

freq = 5; % 2:0.2:20; % 20;

% Grid altitudes
alts = 90:2:400;

% locations in lat lon alt format (deg deg km)
txloc = [36, -58, 0];
rxloc = [36, -76, 0];

% Prop details
OX_mode = 0;
maxdist = 1E2;  % meters from homing
tol = [1e-7 0.01 25];       % ODE solver tolerance and min max stepsizes
refractive_ind = 1;  % at the transmitter location
nhops = 1;

%% Specify time in PHaRLaP format
year = str2double(datestr(time, 'yyyy'));
month = str2double(datestr(time, 'mm'));
day  = str2double(datestr(time, 'dd'));
hour = str2double(datestr(time, 'HH'));
minute = str2double(datestr(time, 'MM'));
UT = [year month day hour minute];


%% Load ionosphere from SAMI3
sami = load_sami(sami_fn);
sami_t = interp_sami(sami, alts, time);

[iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    Bx, By, Bz, geomag_grid_parms] = gen_iono_geomag_grids(sami_t);

iono_pf_grid = sqrt(80.6 * iono_en_grid ./ 1E6);



%% Raytrace
% homed_ray = ...
%     raytrace_itsi(freq, OX_mode, txloc, rxloc, iono_en_grid, iono_en_grid_5, ...
%     collision_freq, iono_grid_parms, Bx, By, Bz, geomag_grid_parms, ...
%     elvarr, azarr, maxdist, tol, 1, 0);


ray = raytrace_3dhome(freq, txloc, rxloc, iono_en_grid, iono_en_grid_5, ...
    collision_freq, iono_grid_parms, Bx, By, Bz, geomag_grid_parms, OX_mode, ...
    tol, refractive_ind, nhops)



%%
plot_rays(txloc, rxloc, ray)



























