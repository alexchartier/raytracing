%% generate a set of vertical ionograms to use for fitting
clear
UT = [2020, 1, 1, 0, 0];
alts = 90:4:600;
lat = 0;
lon = 0;
del_lat = 5;
del_lon = 5;
sc_alt = 580;
foE = 1E6;
hmE = 100;
ymE = 2;
foF1 = 0.1E6;
hmF1 = 220;
ymF1 = 30;
foF2 = 10E6;
hmF2 = 300;
ymF2 = 100;
freqs = 2:0.1:20;

% raytracing stuff
tol = [1e-7 0.01 25];       % ODE solver tolerance and min max stepsizes
nhops = 1;
homing_tol_m = 1500;
OX_mode = 1;
maxdist = 1E5;

%% Generate electron density grid
txloc = [lat, lon, sc_alt];

plasmafreq = chapman(10E6, 100, 2, 0.1E6, 220, 30, foF2, hmF2, 100, alts);
% plot(plasmafreq, alts)
% plot(Ne, alts)
Ne = freq2elec(plasmafreq) / 1E6;

XYZ = igrfmagm(alts * 1E3, lat * ones(size(alts)), lon * ones(size(alts)), 2020 * ones(size(alts)));
Bxi(1, :, :) = XYZ(:, 1)' / 1E9;
Byi(1, :, :) = XYZ(:, 2)' / 1E9;
Bzi(1, :, :) = XYZ(:, 3)' / 1E9;

lats = [lat - del_lat, lat + del_lat];
lons = [lon - del_lon, lon + del_lon];

Ne3(1, :, :) = Ne;
iono_en_grid = repmat(Ne3, length(lats), length(lons), 1);
Bx = repmat(Bxi, length(lats), length(lons), 1);
By = repmat(Byi, length(lats), length(lons), 1);
Bz = repmat(Bzi, length(lats), length(lons), 1);
collision_freq = zeros(size(iono_en_grid));

iono_grid_parms = [...
    min(lats), unique(diff(lats)), length(lats), ...
    min(lons), unique(diff(lons)), length(lons), ...
    min(alts), unique(diff(alts)), length(alts)];

%% Fire rays vertically down
freq = 5;
Ne = interp_3d(lats, lons, alts, iono_en_grid, txloc);

Bx_i = interp_3d(lats, lons, alts, Bx, txloc);
By_i = interp_3d(lats, lons, alts,  By, txloc);
Bz_i = interp_3d(lats, lons, alts,  Bz, txloc);
rsv = calc_rsv(txloc(1), txloc(2), txloc(3), -90, 0, ...
    freq, OX_mode, Ne, Bx_i, By_i, Bz_i);
[~, rays, ~] = ...
    raytrace_3d(txloc(1), txloc(2), txloc(3),  -90, 0, freq, ...
    OX_mode, nhops, tol, iono_en_grid, iono_en_grid, ...
    collision_freq, iono_grid_parms, Bx, By, Bz, ...
    iono_grid_parms, rsv);


%% Home
homed_rays = gen_ionogram(freqs, OX_mode, txloc, txloc, ...
    iono_en_grid, iono_en_grid, collision_freq, ...
    iono_grid_parms, Bx, By, Bz, ...
    maxdist, tol, homing_tol_m);


%% Plot
plot_rays_on_2d_ionosphere(homed_rays, iono_en_grid, lats, lons, alts)


