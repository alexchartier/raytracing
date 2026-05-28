%% OTHR_sim_gitm.m
% Simulate over-the-horizon radar rays using a GITM NetCDF snapshot

gitm_fn = '/Users/chartat1/data/gitm/result_202112_base-001-dene.nc';
out_fn_fmt = '~/data/raytracing/othr_gitm/rays/base/rays_{YYYY-mm-dd_HHMM}.mat';
plt_fn_fmt = fullfile(getenv('HOME'), 'data', 'raytracing', 'othr_gitm', ...
    'plots', 'ba', '{YYYY-mm-dd_HHMM}.png');

% gitm_fn = '/Users/chartat1/data/gitm/result_202112_euv-002-dene.nc';
% out_fn_fmt = '~/data/raytracing/othr_gitm/rays/euv/rays_{YYYY-mm-dd_HHMM}.mat';
% plt_fn_fmt = fullfile(getenv('HOME'), 'data', 'raytracing', 'othr_gitm', ...
%     'plots', 'euv', '{YYYY-mm-dd_HHMM}.png');

gitm_meta = load_gitm_metadata(gitm_fn);
max_snapshots = 306; % increase to process more GITM snapshots
lat_stride = 2;     % keep at 1 for full 1.25 deg resolution
lon_stride = 2;     % keep at 1 for full 2.5 deg resolution
time_indices = 1:min(max_snapshots, numel(gitm_meta.times));
times = gitm_meta.times(time_indices);

ensure_parent_dir(out_fn_fmt);
ensure_parent_dir(plt_fn_fmt);

txloc = [-51.833, -58.983, 0];
freq = 12.2;
OX_mode = 0;
nhops = 3;
tol = [1e-7 0.01 25];

elvarr = 5:5:40;
% azarr = [203.515940660483, 200.140356693153, 196.770909536894, 193.406359932096, 190.045629800606, 186.687760031967, 183.331876773562, 179.9771634617, 176.6228365383, 173.268123226438, 169.912239968033, 166.554370199394, 163.193640067904, 159.829090463106, 156.459643306847, 153.084059339517];
azarr = [203.5159  200.1404  196.7709  193.4064  190.0456  186.6878  183.3319 ...
      179.9772  176.6228  173.2681  169.9122  166.5544  163.1936  159.8291 ...
        156.4596  153.0841];

%% Generate rays
for t = 1:length(times)
    idx = time_indices(t);
    gitm = load_gitm_snapshot(gitm_fn, gitm_meta, idx, lat_stride, lon_stride);
    fprintf('Tracing snapshot %d/%d (%s)\n', t, length(times), datestr(times(t)));
    
    % Generate rays
    rays = gs_raytrace(elvarr, azarr, freq, nhops, OX_mode, ...
        txloc(1), txloc(2), txloc(3), gitm.iono_en_grid, gitm.iono_en_grid, ...
        gitm.collision_freq, gitm.iono_grid_parms, gitm.Bx, gitm.By, gitm.Bz, ...
        gitm.geomag_grid_parms, tol, 0, 0, 0, 0);
    fprintf('Completed raytrace for snapshot %d/%d\n', t, length(times));

    savestruct(filename(out_fn_fmt, times(t)), rays);
    fprintf('Saved to %s\n', filename(out_fn_fmt, times(t)))
end

%% plot
for t= 1:length(times)
    idx = time_indices(t);
    gitm = load_gitm_snapshot(gitm_fn, gitm_meta, idx, lat_stride, lon_stride);

    rays = loadstruct(filename(out_fn_fmt, times(t)));
    txloc(3) = 20;
    % Ionosphere plot
    
    G = flipud(max(gitm.iono_en_grid(2:end-1, :, :), [], 3));

    % G = squeeze(max(gitm.dene, [], 1));
    G(isnan(G)) = 1;
    G = elec2freq(G) * 1E3;
    % simplified jet-like colormap (blue -> cyan -> yellow/green -> orange -> red)
    jet_simple = [ ...
        0.00 0.00 0.35; ... % deep blue
        0.00 0.30 0.80; ... % blue
        0.00 0.65 0.95; ... % cyan
        0.30 0.80 0.50; ... % green-yellow
        0.80 0.88 0.30; ... % yellow
        1.00 0.55 0.12; ... % orange
        0.90 0.10 0.00; ... % red
        ];
    colormap(jet_simple)
    C = colormap;  % Get the figure's colormap.
    L = size(C,1);

    % Scale the matrix to the range of the map.
    Gs = round(interp1(linspace(min(G(:)), max(G(:)), L), 1:L, G));
    Gs(isnan(Gs)) = 1;
    H = reshape(C(Gs(:), :), [size(Gs) 3]);


    space_color = 'k';
    npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
    alpha = 1; % globe transparency level, 1 = opaque, through 0 = invisible

    % Mean spherical earth
    erad    = 6371008.7714; % equatorial radius (meters)
    prad    = 6371008.7714; % polar radius (meters)
    close

    fh = figure('Color', space_color, 'units','normalized','outerposition',[0 0 1 1], ...
        'Position', get(0, 'Screensize')); % fullscreen
    set(gcf, 'WindowState', 'fullscreen');
    earth_example

    hold on;
    % % simple starfield backdrop
    % nstars = 800;
    % star_rad = erad * 10;
    % star_dir = randn(nstars, 3);
    % star_dir = star_dir ./ vecnorm(star_dir, 2, 2);
    % stars = star_dir * star_rad;
    % scatter3(stars(:, 1), stars(:, 2), stars(:, 3), 2, ones(nstars, 3), 'filled', ...
    %     'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 0.0);

    % Turn off the normal axes
    set(gca, 'NextPlot','add', 'Visible','off');
    axis equal;
    axis auto;

    % Set initial view

    axis vis3d;
    view(0, -80); % look down toward the south pole to expose rays


    % Ionospheric structure as semi-transparent isosurfaces
    lats_grid = arange(gitm.iono_grid_parms(1), gitm.iono_grid_parms(2), gitm.iono_grid_parms(3));
    lons_grid = arange(gitm.iono_grid_parms(4), gitm.iono_grid_parms(5), gitm.iono_grid_parms(6));
    alts_grid = arange(gitm.iono_grid_parms(7), gitm.iono_grid_parms(8), gitm.iono_grid_parms(9));
    lon_inc = gitm.iono_grid_parms(5);
    lons_wrap = [lons_grid, lons_grid(1) + lon_inc * numel(lons_grid)];
    lon0_shift = lons_wrap(1) + lon_inc / 2;
    lons_shifted = lon0_shift + (0:numel(lons_wrap)-1) * lon_inc;
    [lat3, lon3, alt3] = ndgrid(lats_grid, lons_shifted, alts_grid);

    iono_grid = gitm.iono_en_grid;
    iono_grid = cat(2, iono_grid, iono_grid(:, 1, :)); % wrap lon to avoid seam
    iono_grid(isnan(iono_grid)) = 0;
    iso_vals = prctile(iono_grid(iono_grid > 0), [70 85 95]);
    iso_vals = iso_vals(iso_vals > 0);

    ax = gca;
    ax.Clipping = "off";

    colormap(C);
    if ~isempty(iso_vals)
        iso_colors = round(linspace(1, size(C, 1), max(numel(iso_vals), 2)));
        iso_alphas = linspace(0.26, 0.46, numel(iso_vals)); % lower densities slightly less transparent
        for iv = 1:numel(iso_vals)
            iso = isosurface(lon3, lat3, alt3, iono_grid, iso_vals(iv));
            sph_iso = [iso.vertices(:, 3) * 1E3 + erad, deg2rad(iso.vertices(:, 2)), deg2rad(iso.vertices(:, 1))];
            iso_cart = sphcart(sph_iso);
            iso.vertices = iso_cart;
            p = patch(iso);
            set(p, 'FaceColor', C(iso_colors(iv), :), 'FaceAlpha', iso_alphas(iv), 'EdgeColor', 'none');
        end
        lighting gouraud
        dt = datetime(times(t), 'ConvertFrom', 'datenum');
        doy = day(dt, 'dayofyear');
        ut_hours = hour(dt) + minute(dt) / 60 + second(dt) / 3600;
        [subsolar_lat, subsolar_lon] = subsolar_point(dt);
        sun_vec = sphcart([erad * 10, deg2rad(subsolar_lat), deg2rad(subsolar_lon)]);
        sun_vec = sun_vec(:)'; % ensure 1x3 row for light position
        sun_vec(1) = sun_vec(1) * 1E30;
        light('Position', sun_vec, 'Style', 'infinite');
    end

    h = colorbar;
    caxis([3e6 8e6]); % 3-8 MHz
    tick_vals = linspace(3e6, 8e6, 6);
    h.Ticks = tick_vals;
    h.TickLabels = compose('%.1f', tick_vals / 1E6);
    h.Color = 'w';
    h.FontSize = 16;
    h.Ruler.Color = 'w';
    h.YColor = 'w';
    h.XColor = 'w';
    h.TickLabelInterpreter = 'tex';
    h.Label.Color = 'w';
    ylabel(h, 'Plasma Frequency (MHz)', 'Color', 'w')
    txt_str = filename('{yyyy-mm-dd HH:MM} UT                                                   ', times(t));
    h.Title.String = txt_str;
    h.Title.Color = 'w';
    h.Title.FontSize = 20;
    h.Title.FontWeight = 'bold';

    hold on
    Re = 6380E3;
    ray_color = [1 1 1];
    highlight_color = 'r';
    highlight_tol_deg = 1;
    % mark transmitter location
    tx_cart = sphcart([txloc(3) * 1E3 + Re, deg2rad(txloc(1)), deg2rad(txloc(2))]);
    plot3(tx_cart(1), tx_cart(2), tx_cart(3), 'r.', 'MarkerSize', 22, 'LineWidth', 2, ...
        'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
    for r = 1:length(rays)

        if ~isempty(rays(r).initial_elev)
            hidx = rays(r).height >= -10;
            sph = [rays(r).height(hidx) * 1E3 + Re; deg2rad(rays(r).lat(hidx)); deg2rad(rays(r).lon(hidx))];
            cart = sphcart(sph');
            plot3(cart(:, 1), cart(:, 2), cart(:, 3), 'Color', ray_color, 'LineWidth', 0.25);

            if isfield(rays, 'wavenorm_B_angle') && numel(rays(r).wavenorm_B_angle) >= numel(hidx)
                ang = rays(r).wavenorm_B_angle(hidx);
                near90 = abs(ang - 90) <= highlight_tol_deg;
                if any(near90)
                    idxNear = find(near90);
                    group_start = idxNear(1);
                    for k = 2:numel(idxNear)
                        if idxNear(k) ~= idxNear(k-1) + 1
                            plot3(cart(group_start:idxNear(k-1), 1), cart(group_start:idxNear(k-1), 2), ...
                                cart(group_start:idxNear(k-1), 3), 'Color', highlight_color, 'LineWidth', 10);
                            group_start = idxNear(k);
                        end
                    end
                    plot3(cart(group_start:idxNear(end), 1), cart(group_start:idxNear(end), 2), ...
                        cart(group_start:idxNear(end), 3), 'Color', highlight_color, 'LineWidth', 2);
                end
            end
        end

        if rays(r).height(end) <=10

            cart = sphcart([rays(r).height(end) * 1E3 + Re, ...
                deg2rad(rays(r).lat(end)), ...
                deg2rad(rays(r).lon(end))]);
            h3 = plot3(cart(1), cart(2), cart(3), 'o', 'markersize', 8, 'markerfacecolor', highlight_color, 'Color', highlight_color);
        end
    end
    hold off

    view_angle = 3.0159548572952004;
    campos = [-38422505.872256, -96019682.10430004, -73115762.6500817];
    camtarget = [309559.73788878, -841730.76691558, -1.0e7]; % aim lower to center south pole vertically
    camup = [0, 0, 1];
    view_angles = [-22.14376826, -72.97123376];

    set(ax, 'CameraPosition', campos, 'CameraTarget', camtarget, ...
        'CameraUpVector', camup, 'CameraViewAngle', view_angle, ...
        'Projection', 'perspective');

    view(ax, view_angles(1), view_angles(2));
    camdolly(0, -0.3, 0)
    export_fig(filename(plt_fn_fmt, times(t)))
    pause(0.1)
    %%
    % clf

end


%% subfunc def
function meta = load_gitm_metadata(gitm_fn)
meta.lat = double(ncread(gitm_fn, 'glat'));
meta.lon = double(ncread(gitm_fn, 'glon'));
meta.alt = double(ncread(gitm_fn, 'galt'));
meta.utc_hours = double(ncread(gitm_fn, 'utc_hrs'));
meta.year = double(ncread(gitm_fn, 'year'));
meta.month = double(ncread(gitm_fn, 'month'));
meta.day = double(ncread(gitm_fn, 'day'));

base_time = datenum(meta.year, meta.month, meta.day);
meta.times = base_time + meta.utc_hours ./ 24;
meta.num_lat = numel(meta.lat);
meta.num_lon = numel(meta.lon);
meta.num_alt = numel(meta.alt);
end


function gitm = load_gitm_snapshot(gitm_fn, meta, idx, lat_stride, lon_stride)
persistent geomag_cache
start = [1, 1, 1, idx];
count = [meta.num_lon, meta.num_lat, meta.num_alt, 1];
dene_slice = permute(ncread(gitm_fn, 'dene', start, count), [3, 2, 1]) ./ 1E6; % m-3 -> cm-3

orig_alt = meta.alt(:);
target_alt = linspace(orig_alt(1), orig_alt(end), meta.num_alt);
if any(abs(diff(orig_alt(1:2)) - diff(orig_alt)) > 1e-6)
    reshaped = reshape(dene_slice, meta.num_alt, []);
    dene_slice = reshape(interp1(orig_alt, reshaped, target_alt, 'linear', 'extrap'), ...
        meta.num_alt, meta.num_lat, meta.num_lon);
    alt_vec = target_alt(:);
else
    alt_vec = orig_alt;
end

dene_full = dene_slice;
lat_idx = 1:lat_stride:meta.num_lat;
lon_idx = 1:lon_stride:meta.num_lon;
dene_slice = dene_slice(:, lat_idx, lon_idx);
lat_vec = meta.lat(lat_idx);
lon_vec = meta.lon(lon_idx);

model.lat = lat_vec(:);
model.lon = lon_vec(:)';
model.alt = alt_vec;
model.time = meta.times(idx);
model.dene = double(dene_slice);

[iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    geomag_grid_parms, grid_axes] = build_iono_grid_struct(model);

if isempty(geomag_cache) || ~isequal(geomag_cache.lat, grid_axes.lat) || ...
        ~isequal(geomag_cache.lon, grid_axes.lon) || ...
        ~isequal(geomag_cache.alt, grid_axes.alt)
    [Bx, By, Bz] = gen_bfield(grid_axes.lat, grid_axes.lon, grid_axes.alt, year(model.time));
    geomag_cache.lat = grid_axes.lat;
    geomag_cache.lon = grid_axes.lon;
    geomag_cache.alt = grid_axes.alt;
    geomag_cache.Bx = Bx;
    geomag_cache.By = By;
    geomag_cache.Bz = Bz;
    geomag_cache.geomag_grid_parms = geomag_grid_parms;
else
    Bx = geomag_cache.Bx;
    By = geomag_cache.By;
    Bz = geomag_cache.Bz;
    geomag_grid_parms = geomag_cache.geomag_grid_parms;
end

gitm.iono_en_grid = iono_en_grid;
gitm.iono_en_grid_5 = iono_en_grid_5;
gitm.collision_freq = collision_freq;
gitm.iono_grid_parms = iono_grid_parms;
gitm.geomag_grid_parms = geomag_grid_parms;
gitm.Bx = Bx;
gitm.By = By;
gitm.Bz = Bz;

gitm.lat = model.lat;
gitm.lon = model.lon;
gitm.lat_full = meta.lat(:);
gitm.lon_full = meta.lon(:)';
gitm.alt = model.alt;
gitm.time = model.time;
gitm.dene = model.dene;
gitm.dene_full = dene_full;
end

function [iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    geomag_grid_parms, grid_axes] = build_iono_grid_struct(model)
hts = model.alt(:)';
lats = model.lat(:)';
lons_in = model.lon(:);

if min(lons_in) <= 0 && max(lons_in) >= 360
    loni = lons_in < 360;
    base_lons = lons_in(loni);
    dene = model.dene(:, :, loni);
    li = base_lons >= 180;
    lons = [base_lons(li) - 360, base_lons(~li)];

    iono_en_grid = permute(dene, [2, 3, 1]);
    iono_en_grid = [iono_en_grid(:, li, :), iono_en_grid(:, ~li, :)];
    iono_en_grid = [iono_en_grid, iono_en_grid(:, 1, :)];
    lons = [lons, lons(1) + 360];
else
    li = lons_in >= 180;
    lons = [lons_in(li) - 360, lons_in(~li)];

    iono_en_grid = permute(model.dene, [2, 3, 1]);
    iono_en_grid = [iono_en_grid(:, li, :), iono_en_grid(:, ~li, :)];
end

iono_en_grid_5 = iono_en_grid;
collision_freq = zeros(size(iono_en_grid));

ht_start = hts(1);
if numel(hts) > 1
    ht_inc = hts(2) - hts(1);
else
    ht_inc = 0;
end
num_ht = numel(hts);

lat_start = lats(1);
if numel(lats) > 1
    lat_inc = lats(2) - lats(1);
else
    lat_inc = 0;
end
num_lat = numel(lats);

lon_start = lons(1);
if numel(lons) > 1
    lon_inc = lons(2) - lons(1);
else
    lon_inc = 0;
end
num_lon = numel(lons);

iono_grid_parms = [lat_start, lat_inc, num_lat, lon_start, lon_inc, num_lon, ...
    ht_start, ht_inc, num_ht];

geom_hts = hts;
if numel(geom_hts) <= 101
    B_ht_inc = ht_inc;
else
    B_ht_inc = (max(geom_hts) - min(geom_hts)) / 50;
    geom_hts = min(geom_hts):B_ht_inc:max(geom_hts);
end
B_lat_inc = lat_inc;
B_lon_inc = lon_inc;

B_ht_start = ht_start;
B_num_ht = numel(geom_hts);
B_lat_start = lat_start;
if B_lat_inc == 0
    B_num_lat = num_lat;
else
    B_num_lat = ceil(num_lat .* lat_inc ./ B_lat_inc);
end
B_lon_start = lon_start;
if B_lon_inc == 0
    B_num_lon = num_lon;
else
    B_num_lon = ceil(num_lon .* lon_inc ./ B_lon_inc);
end

geomag_grid_parms = [B_lat_start, B_lat_inc, B_num_lat, B_lon_start, ...
    B_lon_inc, B_num_lon, B_ht_start, B_ht_inc, B_num_ht];

grid_axes.lat = lats;
grid_axes.lon = lons;
grid_axes.alt = hts;
end


function ensure_parent_dir(fn_fmt)
[folder_path, ~, ~] = fileparts(fn_fmt);
if ~isempty(folder_path) && exist(folder_path, 'dir') == 0
    mkdir(folder_path);
end
end

function [lat, lon] = subsolar_point(dt)
% Calculate subsolar point (lat, lon in degrees, lon east-positive) from datetime
% Uses a simple solar position approximation referenced to J2000 (sufficient for plotting)
% Reference: NOAA/Almanac-style short formulae

% Julian centuries since J2000.0 at 12:00 TT
jd = datenum(dt) + 1721058.5; % convert MATLAB datenum to Julian date
T = (jd - 2451545.0) / 36525.0;

% Days since J2000.0 for mean longitude and anomaly terms
D = jd - 2451545.0;

% Mean longitude and anomaly of the Sun (degrees)
L = 280.460 + 0.9856474 * D;
g = deg2rad(357.528 + 0.9856003 * D);

% Ecliptic longitude and obliquity (radians)
lambda = deg2rad(L + 1.915 * sin(g) + 0.020 * sin(2 * g));
epsilon = deg2rad(23.439 - 0.0000004 * D);

% Declination and right ascension (radians)
delta = asin(sin(epsilon) .* sin(lambda));
alpha = atan2(cos(epsilon) .* sin(lambda), cos(lambda));

% Greenwich Mean Sidereal Time (degrees)
gmst = 280.46061837 + 360.98564736629 * D + 0.000387933 * T.^2 - (T.^3) / 38710000;

% Subsolar longitude (east-positive, degrees) and latitude (degrees)
lon = mod(rad2deg(alpha) - gmst + 180, 360) - 180; % wrap to [-180, 180)
lat = rad2deg(delta);
end
