%% sat_synth_nebula.m
% Simulate two satellite topside sounder mission data

% TODO: local ionosphere grid?

clear

%% file inputs
% base = '/Users/cantrce1/Desktop/Work/Projects/Nebula/';
% out = '/data/osse_9sat_iri_lowsolar/';
out = '/osse_9sat_iri/';
file_fmt = append(base, '/STK_positions/9sat/*.txt');
out_fn_fmt = append(base, out, '/homings/sat%i_sat%i/sim_twosat_rays_%s_%i.mat');
out_iono = append(base, out, '/truth/');

%% runtime parameter definitions
freqs = [2:75]; % 2:0.2:20; % 20;
R12 = 20; %range ~0-200
maxdist = 1E5; % meters from homing
blind_range = 30E3; % two-way (e.g. there-and-back range)
tol = [1e-7, 0.01, 25]; % ODE solver tolerance and min max stepsizes
homing_tol_m = 1500;
load_iono = 0; % ionosphere already generated and save?

% geomagnetic grid increments
B_ht_inc = 4;
B_lat_inc = 2;
B_lon_inc = 4;
I_ht_inc = 4;
I_lat_inc = 2;
I_lon_inc = 4;
pad_lat = 5;
pad_alt = 100;

% Time for IRI
time = datenum(2024, 10, 3, 0, 0, 0);

% times to run (for the orbits)
timeidx = 1:100:1200;

%% read STK satellite position files
files = dir(file_fmt);
for i = 1:size(files, 1)
    fn = join([files(i).folder, '/', files(i).name]);
    fprintf('reading %s\n', fn)
    A = readtable(fn);
    f = join([string(A{:, 1}), A{:, 2}, string(A{:, 3}), string(A{:, 4})]);
    sat{i}.time = datenum(join([string(A{:, 1}), A{:, 2}, string(A{:, 3}), string(A{:, 4})]));
    sat{i}.lat = A{:, 5};
    sat{i}.lon = A{:, 6};
    sat{i}.alt = A{:, 7};
end

assert(length(files) > 0, 'No satellite files found!')

%% define output file format
if ~exist(out_fn_fmt, 'dir')
    mkdir(out_fn_fmt)
end
if ~exist(out_iono, 'dir')
    mkdir(out_iono)
end

%% generate ionosphere
[iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    Bx, By, Bz, geomag_grid_parms, lats, lons, alts] = ...
    gen_ionosphere(sat, time, timeidx, R12, out_iono, B_ht_inc, ...
    B_lat_inc, B_lon_inc, I_ht_inc, ...
    I_lat_inc, I_lon_inc, pad_lat, ...
    pad_alt, load_iono);

%% start raytracing
fprintf('Started a new raytracing expt\n')
j = 0;

%% Loop over satellite pairings
for p1 = 1:size(files, 1) - 1
    for p2 = p1 + 1:size(files, 1)
        fprintf('P1: %i, P2: %i\n', p1, p2)

        %% loop over times
        for ti = 1:length(timeidx)
            t = timeidx(ti);
            fprintf('t: %i\n', t)

            %% Specify transmitter and receiver locs
            txloc = [sat{p1}.lat(t), sat{p1}.lon(t), sat{p1}.alt(t)];
            rxloc = [sat{p2}.lat(t), sat{p2}.lon(t), sat{p2}.alt(t)];

            %% Calculate foF2 for reference
            iono_en_prof = nan(size(alts));
            for i = 1:length(alts)
                iono_en_prof(i) = interp2(lats, lons, ...
                    squeeze(iono_en_grid(:, :, alts == alts(i)))', ...
                    txloc(1), txloc(2));
            end
            fprintf('FoF2: %1.2f\n', sqrt(80.6*max(iono_en_prof)./1E6))

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

                out_fn = sprintf(out_fn_fmt, p1, p2, OX_mode_name, t);

                %% Raytrace
                homed_rays = gen_ionogram(freqs, OX_mode, txloc, rxloc, ...
                    iono_en_grid, iono_en_grid_5, collision_freq, ...
                    iono_grid_parms, Bx, By, Bz, ...
                    maxdist, tol, homing_tol_m);

                %% save
                % clear raytrace_3d
                if size(homed_rays, 2) > 0
                    homed_rays{1}(1).iono_en_grid = iono_en_grid;
                    homed_rays{1}(1).iono_lat = lats;
                    homed_rays{1}(1).iono_lon = lons;
                    homed_rays{1}(1).iono_alt = alts;

                    savestruct(out_fn, homed_rays)
                    fprintf("saved %i homed_rays to %s\n", length(homed_rays), out_fn)

                    %%  construct data array to save as xls
                    for k = 1:size(homed_rays, 2)
                        for l = 1:length(homed_rays{k})
                            j = j + 1;
                            rxsat(j) = p1;
                            rxlat(j) = txloc(1);
                            rxlon(j) = txloc(2);
                            rxalt(j) = txloc(3);
                            txsat(j) = p2;
                            txlat(j) = rxloc(1);
                            txlon(j) = rxloc(2);
                            txalt(j) = rxloc(3);
                            date(j) = string(datestr(time));
                            freq(j) = homed_rays{k}(l).frequency;
                            group_range(j) = homed_rays{k}(l).group_range_to_rx;
                            received_elv(j) = homed_rays{k}(l).initial_elev;
                            received_bearing(j) = homed_rays{k}(l).initial_bearing;
                            geo_dis_to_rx(j) = homed_rays{k}(l).geometric_dist_to_rx;
                            total_absorption(j) = homed_rays{k}(l).total_absorption;
                            perigee_height(j) = homed_rays{k}(l).perigee;
                            refractive_index(j) = homed_rays{k}(l).refractive_index(1);
                            if homed_rays{k}(l).OX_mode == -1
                                mode(j) = "X";
                            elseif homed_rays{k}(l).OX_mode == 1
                                mode(j) = "O";
                            end
                        end
                    end
                end
            end
        end
    end
end

%% save xls file
rxsat = rxsat';
rxlat = rxlat';
rxlon = rxlon';
rxalt = rxalt';
txsat = txsat';
txlat = txlat';
txlon = txlon';
txalt = txalt';
date = date';
freq = freq';
group_range = group_range';
received_elv = received_elv';
received_bearing = received_bearing';
geo_dis_to_rx = geo_dis_to_rx';
total_absorption = total_absorption';
refractive_index = refractive_index';
perigee_height = perigee_height';
mode = mode';
T = table(rxsat, rxlat, rxlon, rxalt, txsat, txlat, txlon, txalt, geo_dis_to_rx, date, freq, group_range, ...
    received_elv, received_bearing, perigee_height, total_absorption, refractive_index, mode);
writetable(T, join([base, out, 'homings_table.xls']))
