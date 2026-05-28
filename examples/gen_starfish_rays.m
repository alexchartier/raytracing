%% gen_starfish_rays.m
% Reproduce the strange starfish/teepee/bird rays
%%%%%%%%%%%%%
% NOTE: There is now a workaround in raytrace_itsi that prevents this
% behavior (modifies txloc lat/lon when a point on a gridline is called). 
% Disable the workaround to allow these rays to occur. 
%%%%%%%%%%%%%
clear

%% inputs
time = datenum(2003, 9, 21, 0, 0, 0);
freqs = 2:0.1:20; % 2:0.2:20; % 20;
rx_tx_dist = 5; % degrees
satlon = -75;
satlon(satlon > 180) = satlon(satlon > 180) - 360;
satalt = 580;

% ionosphere grid
lats = 40:2:70;
lons = -80:5:-70;
alts = 150:3:610;
lons(lons > 180) = lons(lons > 180) - 360;
R12 = 150; %range ~0-200

% elevation and azimuths
elvarr = -80:2:0; %-80:1:0
azarr = 0:5:360; %-10:5:10; %0:5:360

% pharlap parameters
maxdist = 1E5;  % meters from homing
blind_range = 30E3; % two-way (e.g. there-and-back range)
tol = [1e-7 0.01 25];       % ODE solver tolerance and min max stepsizes
homing_tol_m = 1500;
nhops = 1;

txloc = [46 -75 580];
OX_mode = -1;
fprintf('Started a new raytracing expt\n')


%% Generate ionosphere
B_ht_inc = unique(diff(alts));
B_lat_inc = unique(diff(lats));
B_lon_inc = unique(diff(lons));
[iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    Bx, By, Bz, geomag_grid_parms] = gen_iono_geomag_grids(...
    alts, lats, lons, time, R12, B_ht_inc, B_lat_inc, B_lon_inc);

%% Generate ionogram
vert_rays = gen_ionogram(freqs, OX_mode, txloc, txloc, ...
    iono_en_grid, iono_en_grid_5, collision_freq, ...
    iono_grid_parms, geomag_grid_parms, Bx, By, Bz, maxdist, tol, homing_tol_m);



%% Plot
figure
plot_ionogram(vert_rays)
figure
plot_rays_on_2d_ionosphere(vert_rays, iono_en_grid, lats, lons, alts)


%% Print
fprintf('Freq, Elev, Az\n')
for f = 1:length(vert_rays)
    for r = 1:length(vert_rays{f})

        fprintf('%i %i %1.1f %1.4f %1.4f %1.1f\n', f, r, vert_rays{f}(r).frequency, ...
            vert_rays{f}(r).initial_elev, vert_rays{f}(r).initial_bearing, ...
            vert_rays{f}(r).total_absorption)

    end
end

































