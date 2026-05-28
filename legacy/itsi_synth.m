%% itsi_synth.m
% Simulate ITSI mission data

%% inputs

coords_fn = 'iss_coords.txt';
gain_fn = 'realized_gains.txt';
out_fn = 'sim_ionograms.mat';

freqs = 2:0.2:20; % 20;
alt = 150:2:420;
alt_iss = 410;
elvarr = [-90:-80];
azarr = -180:10:180;
R12 = 100;
UT = [2011 9 21 0 0];
OX = -1:1;
OX_names = {'X-mode', 'No B', 'O-mode'};
kp = 3;
maxdist = 3E3;
blind_range = 30E3; % two-way (e.g. there-and-back range)
tol = [1e-7 0.01 25];       % ODE solver tolerance and min max stepsizes

txpower_db = 30;
sigproc_gain = 7;


%% load gain
sim_gain_ascii = asciiread(gain_fn);
gain = str2num(sim_gain_ascii(2:end, :));
txant_gain = interp1(gain(:, 1)/1E6, gain(:, 2), freqs);
rxant_gain = interp1(gain(:, 1)/1E6, gain(:, 3), freqs);


%% load coordinates
crd_txt = asciiread(coords_fn);
times = datenum(crd_txt(:, 1:20));
crd = str2num(crd_txt);
crd = crd(:, 2:3);

times = times(1:10:end);
crd = crd(1:10:end, :);

hour = (times - floor(times)) * 24;
minute = floor((hour - floor(hour)) * 60);
hour = floor(hour);


%% Calculate rays
nanarr = nan(length(times), length(freqs), length(OX));
range = nanarr;
dist = nanarr;
fspl_db = nanarr;
path_loss_db = nanarr;
absorption_db = nanarr;

fof2 = nan(size(times));
hmf2 = nan(size(times));


for t = 1:length(times)
    %% Load ionosphere
    UT(4) = hour(t);
    UT(5) = minute(t);
    lat = round(crd(t, 1));
    lon = round(crd(t, 2));

    lats = lat-2:lat+2;
    lons = lon-2:lon+2;
    lons(lons < -180) = lons(lons < -180) + 360;

    [iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
        Bx, By, Bz, geomag_grid_parms] = ...
        gen_iono_geomag_grids(alt, lats, lons, UT, R12);

    loc = [lat, lon, alt_iss];

    disp(datestr(times(t)))

    %% Get plasma freq and refractive index at the txloc
    iono_pf_grid = sqrt(80.6 * iono_en_grid ./ 1E6);
    plasmafreq = iono_pf_grid(lats == lat, lons == lon, alt == alt_iss);
    assert(length(plasmafreq) == 1, 'Come back and fix this - make a real interpolator')
    prof = iono_pf_grid(lats == lat, lons == lon, :);
    fof2(t) = max(prof);
    hmf2(t) = alt(prof == fof2(t));

    %% Raytrace each frequency that's above the local plasma frequency
    for f = 1:length(freqs)
        fprintf('Freq: %1.1f MHz\n', freqs(f))
        if plasmafreq > freqs(f)
            %disp('local plasma environment > transmit freq - skipping')
            continue
        end

        % Refractive index calculation.
        % Estimate the refractive  index at ray origin. Note this is only the
        % formula for the no mag. field case. Recommend the use of the full
        % Appleton-Hartree equation to calculate the refractive index for the
        % O and X modes
        refractive_ind = sqrt(1 - plasmafreq .^ 2 ./ freqs(f) .^ 2);  % Collisionless, no-B refractive index
        fprintf('Refractive Index: %1.3f \n', refractive_ind)


        %% Raytrace the problem
        for OX_i = 1:length(OX)
            fprintf('OX mode %i\n', OX(OX_i))

            homed_ray = raytrace_itsi(freqs(f), OX(OX_i), loc, loc, ...
                iono_en_grid, iono_en_grid_5, collision_freq, ...
                iono_grid_parms, Bx, By, Bz, geomag_grid_parms, ...
                elvarr, azarr, maxdist, tol, refractive_ind, blind_range, true);

            [range(t, f, OX_i), dist(t, f, OX_i), path_loss_db(t, f, OX_i), ...
                fspl_db(t, f, OX_i), ] = calc_fspl(ray, loc, OX_mode, nhops, tol, ...
                iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
                Bx, By, Bz, geomag_grid_parms, refractive_ind, reflect);

            [di, close_pt, id] = ray_dist(homed_ray, loc, true);

            absorption_db(t, f, OX_i) = homed_ray.absorption(id);



        end
    end
    clear raytrace_3d

end


%% Correct gains
system_gain = txpower_db + txant_gain + rxant_gain + sigproc_gain;
system_gain = repmat(system_gain, [10, 1, 3]);

D = [];
D.freqs = freqs;
D.times = times;
D.path_loss_db = path_loss_db;
D.fspl_db = fspl_db;
D.absorption_db = absorption_db;
D.range = range;
D.system_gain = system_gain;
D.fof2 = fof2;
D.hmf2 = hmf2;
D.lat = crd(:, 1);
D.lon= crd(:, 2);

savestruct(out_fn, D);


%% Reload
D = load(out_fn);
range = D.range;
times = D.times;
absorption_db = D.absorption_db;



%% Ionogram plot
timeidx = 7;
rx_sig_strength_db = D.system_gain + D.path_loss_db - D.absorption_db;
close
hold on
scatter(freqs, squeeze(range(timeidx, :, 1)) ./ 2, [], squeeze(rx_sig_strength_db(timeidx, :, 1)), 'X')
scatter(freqs, squeeze(range(timeidx, :, 2)) ./ 2, [], squeeze(rx_sig_strength_db(timeidx, :, 2)), 'filled', 'd')
scatter(freqs, squeeze(range(timeidx, :, 3)) ./ 2, [], squeeze(rx_sig_strength_db(timeidx, :, 3)), 'filled')

ylim([0, 800])
xlim([0, max(freqs) + 2])
ylabel('virtual depth')
xlabel('Tx Freq (MHz)')
set(gca, 'YDir','reverse')
legend({'X', 'No B', 'O'})
title(sprintf('foF2: %1.1f MHz, hmF2: %1.1f km, S/C alt: %1.1f km', D.fof2(timeidx), D.hmf2(timeidx), alt_iss))
grid on
grid minor
hC = colorbar;
hC.Label.String = 'Received Power (dBm)';
set(gca, 'FontSize', 20, 'CLim', [-120 -80])

%% Time-Frequency-Intensity pcolor plot
% over-plot fof2 line
rx_sig_strength_db = D.system_gain + D.path_loss_db - D.absorption_db;
clf
OX_i = 2;
fti = nan(length(times), length(freqs));
hold on
plot(times + 5/60/24, D.fof2, 'mx-', 'LineWidth', 3)
h = pcolor(times, freqs, rx_sig_strength_db(:, :, OX_i)');
plot(times + 5/60/24, D.fof2, 'mx-', 'LineWidth', 3, 'MarkerSize', 15)

hC = colorbar;
hC.Label.String = 'Received Power (dBm)';
set(gca, 'FontSize', 20, 'CLim', [-120 -80], 'XLim', [times(1), times(end)])
set(h, 'EdgeColor', 'None')

legend({'FoF2', 'Observed'})
title(OX_names{OX_i})

hold off
datetick('keeplimits')
grid on
grid minor







