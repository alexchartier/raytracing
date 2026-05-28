%% gen_vert_ionograms.m
% generate an input set of ionograms

%% Inputs
in_fn_fmt = '/Users/chartat1/data/sami3/2017_tid/sami_mat/{YYYY-mm-dd_HHMM}.mat';
out_fn_fmt = ['/Users/chartat1/data/sami3/gs_ionograms_3d/hiamcm/', ...
    '{YYYY-mm-dd_HHMM}/vert/ionogram_%1.1fN_%1.1fE_%ikm.mat'];
times = datenum(2017, 1, 11, 12, 0, 0):3/24:datenum(2017, 1, 13);
OX_mode = 1;
freqs = 2:0.1:15;

% search space for ionogram generation
txalt = 580; % 500:10:600;
% txlat = 5:10:85;
% txlon = 5:15:355;

txlat = 30:5:60;
txlon = 5:30:360; % :285; % NOTE: don't go too close to the edge (180 line)


% pharlap parameters
maxdist = 1E5;  % meters from homing
tol = [1e-7 0.01 25];       % ODE solver tolerance and min max stepsizes
homing_tol_m = 1500;
fprintf('Started a new raytracing expt\n')


%% Load and reformat SAMI

for t = 1:length(times)
    time = times(t);
sami = loadstruct(filename(in_fn_fmt, time));
sami.lat = sami.lat(2:end-1);
sami.dene = sami.dene(:, 2:end-1, :);
alts = sami.alt;

[iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    Bx, By, Bz, geomag_grid_parms] = gen_grid_parms(sami);
lats = arange(iono_grid_parms(1), iono_grid_parms(2), iono_grid_parms(3));
lons = arange(iono_grid_parms(4), iono_grid_parms(5), iono_grid_parms(6));
alts = arange(iono_grid_parms(7), iono_grid_parms(8), iono_grid_parms(9));


for l1 = 1:length(txlat)
    for l2 = 1:length(txlon)
        for ai = 1:length(txalt)
            %%
            txloc = [txlat(l1), txlon(l2), txalt(ai)];
            homed_rays = gen_ionogram(freqs, OX_mode, txloc, txloc, ...
                iono_en_grid, iono_en_grid_5, collision_freq, ...
                iono_grid_parms, geomag_grid_parms, Bx, By, Bz, ...
                maxdist, tol, homing_tol_m);

            [frq, rg] = calc_ionogram(homed_rays);


            if isempty(frq)
                fprintf('No valid rays - skipping\n')
               continue
            end

            % save
            out_fn = sprintf(filename(out_fn_fmt, time), txloc(1), txloc(2), txloc(3));
            savestruct(out_fn, homed_rays)
            fprintf('Saved to %s\n', out_fn)
        end

    end
end
end




