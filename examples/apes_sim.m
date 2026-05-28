%% sim_apes_data.m
% Simulate the APES mission data for a full orbit or more

% TODO: 

%% Set inputs
in_fn_fmt = 'data/sami/{YYYY-mm-dd_HHMM}.mat';
out_fn_fmt = 'data/apes_osse/sim_twosat_rays_%s_%i_%1.1f.mat';
time = datenum(2015, 3, 21, 8, 0, 0);

freqs = reshape(linspace(sqrt(2), sqrt(20), 100).^2, 10, 10);
freqs = freqs([1, 5, 8, 3, 10, 2, 7, 4, 9, 6], :)';
freqs = freqs(:);
OX_modes = [-1, 1]; 

satlat = -90; 
satlon = -42.5;
satalt = 800;
sep = 5; 

elvarr = -70:2:0;
azarr = -10:5:10;

maxdist = 1E5;  % meters from homing
tol = [1e-7 0.01 25];       % ODE solver tolerance and min max stepsizes
homing_tol_m = 750; 

fprintf('Started a new raytracing expt\n')

integration_time = 0.1; 
lat_step_deg = 7.8 * integration_time / 120;   % how far the satellite moves in one integration


%% Load ionosphere
% NOTE: not sure there's any point in reducing the grid sizing in terms of
% speed

fprintf('Loading ionosphere\n')

sami = loadstruct(filename(in_fn_fmt, time));
[iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    Bx, By, Bz, geomag_grid_parms] = gen_grid_parms(sami);
lats = arange(iono_grid_parms(1), iono_grid_parms(2), iono_grid_parms(3));
lons = arange(iono_grid_parms(4), iono_grid_parms(5), iono_grid_parms(6));
alts = arange(iono_grid_parms(7), iono_grid_parms(8), iono_grid_parms(9));


%% loops 
pole_passes = 0;
ct = 55836;
fi = 1;

while pole_passes < 2
    
    %% Specify transmitter and receiver locs
    txloc = [satlat, satlon, satalt];
    rxloc = [satlat + sep, satlon, satalt];

    %% Calculate foF2 for reference
    iono_en_prof = nan(size(alts));
    for i = 1:length(alts)
        iono_en_prof(i) = interp2(lats, lons, ...
            squeeze(iono_en_grid(:, :, alts == alts(i)))',...
            txloc(1), txloc(2));
    end
    fprintf('FoF2: %1.2f\n', sqrt(80.6 * max(iono_en_prof)./ 1E6))

    %% Loop over O/X mode
    for oxi = 1:length(OX_modes)
        OX_mode = OX_modes(oxi);
        switch OX_mode
            case -1
                OX_mode_name = 'X';
            case 1
                OX_mode_name = 'O';

            case 0
                OX_mode_name = 'No_B';
        end

        out_fn = sprintf(out_fn_fmt, OX_mode_name, ct, freqs(fi));


        %% perform raytracing
        ray = ...
            raytrace_itsi(freqs(fi), OX_mode, txloc, rxloc, ...
            iono_en_grid, iono_en_grid_5, collision_freq, ...
            iono_grid_parms, Bx, By, Bz, geomag_grid_parms, ...
            elvarr, azarr, maxdist, tol, 0, homing_tol_m);

        %% save out
        savestruct(out_fn, ray)

    end


    %% increment and check for pole passing
    satlat = satlat + lat_step_deg;

    if abs(satlat) >= 90
        pole_passes = pole_passes + 1;
        lat_step_deg = - lat_step_deg;
        sep = - sep;
        satlat = satlat - (satlat - 90);
        satlon = satlon + 180;
        azarr = azarr + 180;
        satlon(satlon > 180) = satlon(satlon > 180) - 360;

    end
    ct = ct + 1;

    if fi < length(freqs)
        fi = fi + 1;
    else 
        fi = 1;
    end
end


%% plotting

% plot(freqs, '.', 'markersize', 20); grid on; grid minor














