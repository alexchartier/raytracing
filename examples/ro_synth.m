% time = datenum(2017, 1, 10, 18, 0, 0);

%NOTE: need to comment out reflection testing
mod_fn = '/Users/chartat1/data/nebula/model_input/iri_2000km/2012-01-01_0000.mat';

freq = 1500;
txloc_1 = [68.7685866015248, -102.56293471036047, 1983.65184792664];
txloc_2 = [19.6254,    -8.7478, 2000];

rxloc = [71.64, -56.91, 573.685];

OX_mode = -1;
% rxloc = [50.82, -41.18, 567.032];
% txloc = [44.07, -170.57, 21551.0789];


maxdist = 1E7;  % meters from homing
blind_range = 30E3; % two-way (e.g. there-and-back range)
tol = [1e-7 0.01 5];       % ODE solver tolerance and min max stepsizes
homing_tol_m = 1500;

%%
model = loadstruct(mod_fn);

[gc_az, gc_el, sep] = geodetic2aer( ...
    rxloc(1), rxloc(2), rxloc(3), ...
    txloc_1(1), txloc_1(2), txloc_1(3), ...
    wgs84Ellipsoid("km"));

[gc_az2, gc_el2, sep] = geodetic2aer( ...
    rxloc(1), rxloc(2), rxloc(3), ...
    txloc_2(1), txloc_2(2), txloc_2(3), ...
    wgs84Ellipsoid("km"));

%%

r1 = raytrace_itsi(freq, OX_mode, txloc_1, rxloc, ...
    model.iono_en_grid, model.iono_en_grid, model.collision_freq, ...
    model.iono_grid_parms, model.Bx, model.By, model.Bz, ...
    model.geomag_grid_parms, ...
    gc_el, gc_az, maxdist, tol, blind_range, homing_tol_m);

r2 = raytrace_itsi(freq, OX_mode, txloc_2, rxloc, ...
    model.iono_en_grid, model.iono_en_grid, model.collision_freq, ...
    model.iono_grid_parms, model.Bx, model.By, model.Bz, ...
    model.geomag_grid_parms, ...
    gc_el2, gc_az2, maxdist, tol, blind_range, homing_tol_m);

%%
earth_example;
plot_rays(r1, txloc_1, rxloc)
plot_rays(r2, txloc_2, rxloc)


