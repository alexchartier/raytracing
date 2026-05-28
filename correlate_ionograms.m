%% Find the closest match between the observed ionogram and the input set
% See also gen_obs_ionograms.m, preproc_ionograms.m, calc_tdoa.m

% TODO:
% - Expand the global-search sets for vertical and oblique
% - Fix the bottomside (Chapman?)

clear
%% Set inputs
time = datenum(2017, 1, 11, 17, 0, 0);

vert = false;
if vert
    iono_set_fn = '/Users/chartat1/data/sami3/gs_ionogram_set_vert.mat';
    obs_fn_fmt = ['/Users/chartat1/data/sami3/2017_tid/recon_rays/O_mode/', ...
        'pass_{yyyy-mm-dd_HHMM}_vert/{yyyy-mm-dd_HHMM}_%i.mat'];
    out_fn_fmt = ['/Users/chartat1/data/sami3/2017_tid/recon_iono/', ...
        'vert_recon_{YYYY-mm-dd_HHMM}_%1.0fN_%1.0fE.mat'];
else
    iono_set_fn = '/Users/chartat1/data/sami3/gs_ionogram_set_oblique_600.mat';
    obs_fn_fmt = ['/Users/chartat1/data/sami3/2017_tid/recon_rays/O_mode/', ...
        'pass_{yyyy-mm-dd_HHMM}_oblique/{yyyy-mm-dd_HHMM}_%i.mat'];
    out_fn_fmt = ['/Users/chartat1/data/sami3/2017_tid/recon_iono/', ...
        'oblique_recon_{YYYY-mm-dd_HHMM}_%1.0fN_%1.0fE.mat'];
    alongtrack_spacing_km = 600;
end

mod_fn_fmt = '/Users/chartat1/data/sami3/2017_tid/sami_mat/{YYYY-mm-dd_HHMM}.mat';
R12 = 30;
sl_cutoff_pct = 10;  % What proportion of ionogram to cut off (within X percent of straight line)
WGS84 = wgs84Ellipsoid;

alts_out = 90:600;


%% Load model inputs
iono_set = loadstruct(iono_set_fn);
sami = loadstruct(filename(mod_fn_fmt, time));

%% Loop through observation set
nobs = length(dir(filename(fileparts(obs_fn_fmt), time))) - 2;

lats_out = zeros(nobs, 1);
lons_out = zeros(nobs, 1);
alts_out = sami.alt;
hmax_out = zeros(nobs, 1);
nemax_out = zeros(nobs, 1);
H_out = zeros(nobs, 1);
ne_local = zeros(nobs, 1);
truth_profs = zeros(length(sami.alt), nobs);
uncal_profs = zeros(length(sami.alt), nobs);
obs_profs = zeros(length(sami.alt), nobs);

assert(length(nobs) > 0, "need observations to reconstruct")

for i = 1:nobs

    %% Load each observation
    obs_rays = loadstruct(sprintf(filename(obs_fn_fmt, time), i));
    txloc = obs_rays{1}(1).txloc;
    obs_iono = get_ionogram(obs_rays, sl_cutoff_pct);
    obs_iono.mod_prof = interp_sami(sami, obs_iono.lat_ref, obs_iono.lon_ref);
    obs_frange = range(obs_iono.freqs);
    truth_profs(:, i) = obs_iono.mod_prof;

    % Store other stuff
    ne_local(i) = obs_iono.N_local;
    lats_out(i) = obs_iono.lat_ref;
    lons_out(i) = obs_iono.lon_ref;


    %% Go through and calculate range/frequency corrections for the 'close' ones
    max_R = 0;
    overall_best_score = 1E6; % Ideal score is zero
    clear best_fi

    %% Check for good obs
    if isempty(obs_frange)
        fprintf(" Obs %i invalid - skipping\n", i)
        continue
    end

    for fi = 1:length(iono_set)
        %% Only look at the closest matches
        % fprintf('obs_frange: %1.1f mod_frange: %1.1f\n', obs_frange, range(iono_set{fi}.freqs))

        % Skip if frequency extents are not similar
        % TODO: tighten up later, when DB is larger
        if abs(obs_frange - range(iono_set{fi}.freqs)) > 1.5
            continue
        end

        % Skip if required adjustment is >1 MHz
        if abs(max(obs_iono.freqs) - max(iono_set{fi}.freqs)) > 1
            continue
        end

        %% Update the iono_set frequency range to match observed peak density
        fdiff_max = max(obs_iono.freqs) - max(iono_set{fi}.freqs);
        iono_set{fi}.fnew = iono_set{fi}.freqs + fdiff_max;
        iono_set{fi}.fdiff = fdiff_max;

        %% Find the best range adjustment
        best_score = 1E6; % Ideal score is zero
        % Trial range adjustments between -100:100

        for rg_diff = -100:10:100

            %% calculate score for each range correction
            score = 0;
            for fi_mod = 1:length(iono_set{fi}.fnew)
                fi_obs = ismember(...
                    round(obs_iono.freqs * 10), ...
                    round(iono_set{fi}.fnew(fi_mod) * 10));
                if sum(fi_obs) > 0
                    score = score + min(abs(obs_iono.ranges(fi_obs) - ...
                        (iono_set{fi}.ranges(fi_mod) + rg_diff)));
                end
            end

            %% penalize non-matching frequencies
            score = score + 10 * length(...
                setxor(round(iono_set{fi}.fnew * 10), round(obs_iono.freqs * 10)));

            if score < best_score
                best_score = score;
                iono_set{fi}.score = score;
                iono_set{fi}.rnew = iono_set{fi}.ranges + rg_diff;
                iono_set{fi}.rg_diff = rg_diff;
            end
        end

        if best_score < overall_best_score
            best_fi = fi;
            overall_best_score = best_score;
        end
    end

    %% Pull out the parameters and store the winning profile
    mod_iono = iono_set{best_fi};
    uncal_profs(:, i) = mod_iono.mod_prof;

    % Figure out nemax from parameter fitting
    coeffs = iono_set{1}.mdl_fmax_sqrtnmax.Coefficients.Estimate;
    c = coeffs(1);
    m = coeffs(2);
    nemax_out(i) = (m * max(obs_iono.freqs) + c)^2;

    % Range/vht adjustment (may not work that well)
    coeffs = iono_set{1}.mdl_grouprg_peakht.Coefficients.Estimate;
    m = coeffs(2);
    alt_adj = m * mod_iono.rg_diff;
    mod_alts = mod_iono.prof_alts - alt_adj;
    hmax_out(i) = mod_alts(mod_iono.mod_prof == max(mod_iono.mod_prof));
end


%% Fill in gaps
gi = nemax_out == 0;
nemax_out(gi) = interp1(lats_out(~gi), nemax_out(~gi), lats_out(gi), 'linear', 'extrap');
hmax_out(gi) = interp1(lats_out(~gi), hmax_out(~gi), lats_out(gi), "linear", 'extrap');


%% Calculate params
% Store truth vals for posterity
nemax_truth = max(truth_profs, [], 1);
hmax_truth = zeros(size(lats_out));
h_bs_truth = zeros(size(lats_out));
for i = 1:length(lats_out)
    [~, maxi] = max(truth_profs(:, i));
    hmax_truth(i) = sami.alt(maxi);
    H_bs_truth(i) = find_scaleheight(truth_profs(:, i)', sami.alt, false);

end

sc_spacing = dist_txrx(obs_iono.txloc, obs_iono.rxloc);
meas_spacing = dist_txrx([lats_out(1), lons_out(1), 300], [lats_out(2), lons_out(2), 300]);

% Smoothing and calibration
nemax_out_cal = smooth(nemax_out, 15);
nemax_grad = gradient(nemax_out_cal, meas_spacing);
hmax_out_cal = smooth(hmax_out, 15); %  .* abs(1 + (sc_spacing / 600 * nemax_grad / 4) / 100);


%% Generate profiles from the parameters
sc_alt = obs_iono.txloc(3);
for i = 1:length(lats_out)
    obs_profs(:, i) = gen_iri_profile(time, nemax_out_cal(i), hmax_out_cal(i), ...
        ne_local(i), sc_alt, lats_out(i), lons_out(i), alts_out, R12);
end


%% Map to a regular grid and save
out.time = time;
out.alt = alts_out;
out.obs_lat = lats_out;
out.obs_lon = lons_out;
out.profs = obs_profs;
out.truth_profs = truth_profs;
out.uncal_profs = uncal_profs;

% map out into 3D (constant in latitude)
grid_lats = min(lats_out):max(lats_out);
grid_lons = min(lons_out) -10:max(lons_out) + 10;
grid_dene = zeros(length(mod_alts), length(grid_lats), length(grid_lons));
for alti = 1:length(mod_alts)
    dene_i1 = interp1(lats_out, obs_profs(alti, :), grid_lats);
    for loni = 1:length(grid_lons)
        grid_dene(alti, :, loni) = dene_i1;
    end
end

out.dene = grid_dene;
out.lat = grid_lats;
out.lon = grid_lons;

% add the B-field
[iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    Bx, By, Bz, geomag_grid_parms] = gen_grid_parms(out);
out.iono_en_grid = iono_en_grid;
out.collision_freq = collision_freq;
out.iono_grid_parms = iono_grid_parms;
out.geomag_grid_parms = geomag_grid_parms;
out.Bx = Bx;
out.By = By;
out.Bz = Bz;

% save
out_fn = filename(sprintf(out_fn_fmt, mean(out.lat), mean(out.lon)), time);
savestruct(out_fn, out)
fprintf('Saved to %s\n', out_fn)


%% Calibrated params plot
close
subplot(2, 1, 1)
hold on
plot(lats_out, nemax_truth, '-k', 'LineWidth', 2)
plot(lats_out, nemax_out, '-r', 'LineWidth', 2)
plot(lats_out, nemax_out_cal, '--r', 'LineWidth', 2)
legend({'Truth', 'Observed', 'Observed (smoothed)'})


ylabel('NmF2 (el. cm^{-3}')
ylim([0, 4E5])
grid on
grid minor

subplot(2, 1, 2)
hold on
plot(lats_out, hmax_truth, '-k', 'LineWidth', 2)
plot(lats_out, hmax_out, '-r', 'LineWidth', 2)
plot(lats_out, hmax_out_cal, '--r', 'LineWidth', 2)
legend({'Truth', 'Observed (uncal)', 'Smoothed & calibrated'})
ylim([0, 400])
xlabel('Lat (°)')
ylabel('hmF2 (km)')
grid on
grid minor

figure
hold on
i = 8;
plot(truth_profs(:, i), alts_out, '-k')
plot(obs_profs(:, i), alts_out, '-r')
title(sprintf('%s UT %1.1° N', datestr(time), lats_out(i)))
xlabel("Electron Density (el. cm-3")
ylabel("Alt (km)")
grid on
grid minor
hold off
