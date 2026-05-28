%% Find the closest match between the observed ionogram and the input set
% See also gen_obs_ionograms.m, preproc_ionograms.m, calc_tdoa.m

% TODO: 
% - Expand the global-search sets for vertical and oblique
% - Fix the bottomside (Chapman?)
% - Map to the across-track direction properly

clear
%% Set inputs
vert = false;
if vert
    iono_set_fn = '/Users/chartat1/data/sami3/gs_ionogram_set_vert.mat';
    obs_fn_fmt = ['/Users/chartat1/data/sami3/2017_tid/recon_rays/O_mode/', ...
        'vert_{yyyy-mm-dd_HHMM}_%1.1fN_%1.1fE_%ikm.mat'];
    out_fn_fmt = ['/Users/chartat1/data/sami3/2017_tid/recon_iono/', ...
        'vert_recon_{YYYY-mm-dd_HHMM}_%1.1fN_%1.1fE.mat'];
else
    iono_set_fn = '/Users/chartat1/data/sami3/gs_ionogram_set_oblique_600.mat';
    obs_fn_fmt = ['/Users/chartat1/data/sami3/2017_tid/recon_rays/', ...
    'O_mode/oblique_{YYYY-mm-dd_HHMM}_%1.1fN_%1.1fE_%1.1fN_%1.1fE_%ikm.mat'];
    out_fn_fmt = ['/Users/chartat1/data/sami3/2017_tid/recon_iono/', ...
        'oblique_recon_{YYYY-mm-dd_HHMM}_%1.0fN_%1.0fE.mat'];
    alongtrack_spacing_km = 600;
end

mod_fn_fmt = '/Users/chartat1/data/sami3/2017_tid/sami_mat/{YYYY-mm-dd_HHMM}.mat';
time = datenum(2017, 1, 10, 18, 0, 0);
lats = 30:55;
lon = -77.5;
sc_alt = 580;

WGS84 = wgs84Ellipsoid; 

%% Load model inputs
iono_set = loadstruct(iono_set_fn);
sami = loadstruct(filename(mod_fn_fmt, time));

%% Loop through observation set
profs = zeros(length(sami.alt), length(lats)-1 );
obs_profs = zeros(length(sami.alt), length(lats)- 1);

lats_out = zeros(size(lats) - 1);
lons_out = zeros(size(lats)- 1);
fof2 = zeros(size(lats)-1);
n_local = zeros(size(lats) - 1);
elv = zeros(size(lats)-1);
peak_vht_diff = zeros(size(lats)-1);
for i = 1:length(lats) - 1

    %% Load each observation
    txloc = [lats(i), lon, sc_alt];
    if vert
        obs_rays = loadstruct(sprintf(filename(obs_fn_fmt, time), txloc(1), txloc(2), txloc(3)));
    else
        refloc = [lats(i + 1), lon, sc_alt];
        rxloc = calc_rx_loc_for_ois(txloc, refloc, alongtrack_spacing_km);
        obs_rays = loadstruct(sprintf(filename(obs_fn_fmt, time), ...
                    txloc(1), txloc(2), rxloc(1), rxloc(2), txloc(3)));
    end
    
    obs_iono = cleanup_ionogram(obs_rays);
    lats(i) = obs_iono.txloc(1);
    fof2(i) = max(obs_iono.freqs);

    %% loop through ionogram set
    min_chi2 = 1E9;
    idx = [];
    for fi = 1:length(iono_set)

        %% Check the frequency range is within 1 MHz
        if abs(range(iono_set{fi}.freqs) - range(obs_iono.freqs)) > 1
            continue
        end

        %% Match the peak frequencies
        iono_set{fi}.fdiff = max(obs_iono.freqs) - max(iono_set{fi}.freqs);
        iono_set{fi}.fnew = iono_set{fi}.freqs + iono_set{fi}.fdiff;

        %% Determine agreement
        % get the matching frequency indices
        fidx_m = ismember(round(iono_set{fi}.fnew * 10), round(obs_iono.freqs * 10));
        fidx_o = ismember(round(obs_iono.freqs * 10), round(iono_set{fi}.fnew * 10));

        vht_diffs = iono_set{fi}.vht(fidx_m) - obs_iono.vht(fidx_o);
        chi2 = sum(vht_diffs .^ 2)  / length(vht_diffs);
        if chi2 < min_chi2
            min_chi2 = chi2;
            idx = fi;
        end

    end
    mod_iono = iono_set{idx};
    peak_vht_diff(fi) = max(obs_iono.vht) - max(mod_iono.vht);

    %% figure out where the measurement is located
    % save out 'actual' locations, using observed azimuth plus distance
    % between transmitter and reflected location from selected model ionogram
    mod_refloc = [mean(mod_iono.lat_ref(end-5:end)), ...
        mean(mod_iono.lon_ref(end-5:end)), mod_iono.txloc(3)];
    
    arclen_m = distance(mod_iono.txloc(1), mod_iono.txloc(2),...
        mod_refloc(1), mod_refloc(2), WGS84, 'degrees');
    az = mean(obs_iono.az(end-5:end));
    [lats_out(i), lons_out(i)] = reckon(txloc(1), txloc(2), arclen_m, az, WGS84, 'degrees');

    %% Pull out a profile and rescale it
    sami = loadstruct(filename(mod_fn_fmt, time));

    mod_prof = interp_sami(sami, mod_refloc);
    fmax_obs = max(obs_iono.freqs);
    fmax_mod = max(mod_iono.freqs);
    fof2_fac = predict(iono_set{1}.mdl_fmax_fof2, fmax_obs) - ...
        predict(iono_set{1}.mdl_fmax_fof2, fmax_mod);

    mod_prof_rescale = freq2elec((elec2freq(mod_prof) / 1E3) + fof2_fac) * 1E6;
    profs(:, i) = mod_prof_rescale;

    elv(i) = mean(obs_iono.el(end-5:end));

    n_local(i) = obs_iono.N_local;
end


%% Map to a regular grid and save
out.time = time;
out.alt = sami.alt;
out.obs_lat = lats_out;
out.obs_lon = lons_out;
out.profs = profs; 

lati = ceil(min(lats_out)):floor(max(lats_out) - 2);
dene = zeros(length(out.alt), length(lati));
for i = 1:length(out.alt)
    dene(i, :) = interp1(out.obs_lat, out.profs(i,:), lati);
end

% scale to match in situ density observations
alts = squeeze(repmat(out.alt, [1, size(lati)]));
hmf2 = alts(dene == max(dene));
for i = 1:length(lati)
    topside_idx = out.alt >= hmf2(i);
    scale_fac = n_local(round(lats) == lati(i)) / dene(out.alt == sc_alt, i);
    scaling = zeros(sum(topside_idx), 1);
    peak_to_sc_idx = out.alt >= hmf2(i) & out.alt <= sc_alt;
    scaling(1:sum(peak_to_sc_idx)) = linspace(1, scale_fac, sum(peak_to_sc_idx));
    scaling(sum(peak_to_sc_idx)+1:end) = interp1(...
        out.alt(peak_to_sc_idx), scaling(1:sum(peak_to_sc_idx)), ...
        out.alt(out.alt > sc_alt), 'linear', 'extrap');
    dene(topside_idx, i) = dene(topside_idx, i) .* scaling;    
end


% Replace bottomside with Chapman
dene0 = dene; 
H = zeros(size(hmf2));
for i = 1:size(dene, 2)
    H(i) = find_scaleheight(dene(:, i), out.alt);
end
H(end) = H(end - 1);
H = smooth(H, 3);
for i = 1:size(dene, 2)
dene(out.alt < hmf2(i), i) = ...
    calc_chapman_prof(out.alt(out.alt < hmf2(i)), max(dene(:, i)), hmf2(i), H(i), 1);
end
% map out into 3D
out.dene = repmat(dene, [1, 1, 3]);
out.lat = lati;
out.lon = [mean(lons_out) - 2, mean(lons_out), mean(lons_out) + 2];


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

