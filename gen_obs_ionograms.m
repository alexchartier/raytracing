%% gen_obs_ionograms.m
% Generate a simulated 'observed' set of oblique ionograms in flyby of
% ground transmitter site

%TODO: Update to use correct coords


clear
%% inputs
sami_fn_fmt = '/Users/chartat1/data/sami3/2017_tid/sami_mat/{YYYY-mm-dd_HHMM}.mat';
% sc_fn_fmt = ['~/data/nebula/STK_positions/Nebula_Constellation_LLA_AllTimesteps/', ...
%     'Nebula_Constellation_LLA_Time%i.csv'];

sc_crd_fn1 = '~/data/nebula/STK_positions/tdoa_sats//CubeSat1_Fixed_Position_Velocity.txt';
sc_crd_fn2 = '~/data/nebula/STK_positions/tdoa_sats//CubeSat2_Fixed_Position_Velocity.txt';

out_fn_fmt_oblique = ['/Users/chartat1/data/sami3/2017_tid/recon_rays/O_mode/', ...
    'pass_{yyyy-mm-dd_HHMM}_oblique/{yyyy-mm-dd_HHMM}_%i.mat'];
out_fn_fmt_vert = ['/Users/chartat1/data/sami3/2017_tid/recon_rays/O_mode/', ...
    'pass_{yyyy-mm-dd_HHMM}_vert/{yyyy-mm-dd_HHMM}_%i.mat'];

txloc_gd = [38, -77.5, 0];
time = datenum(2017, 1, 11, 23, 0, 0);
OX_mode = 1;
freqs = 2:0.1:20;

alongtrack_spacing_km = 600;  % spacing between spacecraft

% pharlap parameters
maxdist = 1E5;  % meters from homing
tol = [1e-7 0.01 25];       % ODE solver tolerance and min max stepsizes
homing_tol_m = 1500;
fprintf('Started a new raytracing expt\n')
gen_vert = true;

wgs84 = wgs84Ellipsoid;

%% Read coordinate data
data_sc1 = readtable(sc_crd_fn1);
data_sc2 = readtable(sc_crd_fn2);
sc_idx = 500:10:700;

sep = 54;
XYZ_sc1 = [data_sc1.Var5(sc_idx), data_sc1.Var6(sc_idx), data_sc1.Var7(sc_idx)] * 1E3;
XYZ_sc2 = [data_sc2.Var5(sc_idx + sep), data_sc2.Var6(sc_idx + sep), data_sc2.Var7(sc_idx + sep)] * 1E3;


%% Create transmitter and receiver LLA locations
SphV = cartsph(XYZ_sc1);
sc_rad = mean(SphV(:, 1));
alt = sc_rad /1E3 - 6371;
txlats = rad2deg(SphV(:, 2));
txlons = rad2deg(SphV(:, 3));

SphV = cartsph(XYZ_sc2);
rxlats = rad2deg(SphV(:, 2));
rxlons = rad2deg(SphV(:, 3));


%% Load and reformat SAMI 
sami = loadstruct(filename(sami_fn_fmt, time));
lons = sami.lon;
lons(lons > 180) = lons(lons > 180) - 360;

%% Loop through tx/rx locations
for l1 = 1:length(txlats)
    %% Generate oblique ionograms and save
    txloc = [txlats(l1), txlons(l1), alt];
    rxloc = [rxlats(l1), rxlons(l1), alt];
    
    homed_rays = gen_ionogram(freqs, OX_mode, txloc, rxloc, ...
        sami.iono_en_grid, sami.iono_en_grid, sami.collision_freq, ...
        sami.iono_grid_parms, sami.geomag_grid_parms, sami.Bx, sami.By, sami.Bz, ...
        maxdist, tol, homing_tol_m);

    [frq, rg] = calc_ionogram(homed_rays);

    if isempty(frq)
        fprintf('No valid rays - skipping\n')
        continue
    end
    % save
    out_fn = sprintf(filename(out_fn_fmt_oblique, time), l1);
    savestruct(out_fn, homed_rays) 

    %% Generate vertical ionograms and save
    if gen_vert
        txloc = [txlats(l1), txlons(l1), alt];
        rxloc = txloc;

        homed_rays = gen_ionogram(freqs, OX_mode, txloc, rxloc, ...
            sami.iono_en_grid, sami.iono_en_grid, sami.collision_freq, ...
            sami.iono_grid_parms, sami.geomag_grid_parms, sami.Bx, sami.By, sami.Bz, ...
            maxdist, tol, homing_tol_m);

        [frq, rg] = calc_ionogram(homed_rays);

        if isempty(frq)
            fprintf('No valid rays - skipping\n')
            continue
        end
        % save
        out_fn_vert = sprintf(filename(out_fn_fmt_vert, time), l1);
        savestruct(out_fn_vert, homed_rays)
        fprintf('Saved to %s\n', out_fn_vert)
    end
end



































