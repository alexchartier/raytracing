%% twosat_synth.m
% Simulate two satellite topside sounder mission data

clear

%% inputs
out_fn_fmt = 'data/iri_osse_el80-0_az0-360/sim_twosat_rays_%s_%i.mat';

%% runtime parameter definitions
time = datenum(2003, 9, 21, 0, 0, 0);
freqs = 3:0.1:20; % 2:0.2:20; % 20;
satlats = 40:65;  % will move through these
rx_tx_dist = 5; % degrees
satlon = -75;
satlon(satlon > 180) = satlon(satlon > 180) - 360;
satalt = 600;

% ionosphere grid
lats = 40:2:70;
lons = -80:2:-70;
alts = 150:2:610;
lons(lons > 180) = lons(lons > 180) - 360;
R12 = 150; %range ~0-200
B_ht_inc = 10;
B_lat_inc = 10;
B_lon_inc = 10;

% pharlap parameters
maxdist = 1E5;  % meters from homing
tol = [1e-7 0.01 25];       % ODE solver tolerance and min max stepsizes
homing_tol_m = 1500;
fprintf('Started a new raytracing expt\n')

%% Generate ionosphere
[iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    Bx, By, Bz, geomag_grid_parms] = gen_iono_geomag_grids(...
    alts, lats, lons, time, R12, B_ht_inc, B_lat_inc, B_lon_inc);


%% loop over latitude
for lati = 1:range(satlats) - rx_tx_dist
    %% Specify transmitter and receiver locs
    txloc = [satlats(lati), satlon, satalt];
    rxloc = [satlats(lati + rx_tx_dist), satlon, satalt];

    %% Calculate foF2 for reference
    iono_en_prof = nan(size(alts));
    for i = 1:length(alts)
        iono_en_prof(i) = interp2(lats, lons, ...
            squeeze(iono_en_grid(:, :, alts == alts(i)))',...
            txloc(1), txloc(2));
    end
    fprintf('FoF2: %1.2f\n', sqrt(80.6 * max(iono_en_prof)./ 1E6))

    %% Loop over O/X mode
    for OX_mode = -1:2:1

        switch OX_mode
            case -1
                OX_mode_name = 'X';
            case 1
                OX_mode_name = 'O';

            case 0
                OX_mode_name = 'No_B';
        end

        out_fn = sprintf(out_fn_fmt, OX_mode_name, lati);

        %% Raytrace

        homed_rays = gen_ionogram(freqs, OX_mode, txloc, rxloc, ...
            iono_en_grid, iono_en_grid_5, collision_freq, ...
            iono_grid_parms, Bx, By, Bz, ...
            maxdist, tol, homing_tol_m);

    end

    %% save
    clear raytrace_3d
    homed_rays{1}(1).iono_en_grid = iono_en_grid;
    homed_rays{1}(1).iono_lat = lats;
    homed_rays{1}(1).iono_lon = lons;
    homed_rays{1}(1).iono_alt = alts;

    savestruct(out_fn, homed_rays)
    fprintf("saved %i homed_rays to %s\n", length(homed_rays), out_fn)

end


