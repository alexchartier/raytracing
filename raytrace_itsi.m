function homed_ray = ...
    raytrace_itsi(freq, OX_mode, txloc, rxloc, iono_en_grid, iono_en_grid_5, ...
    collision_freq, iono_grid_parms, Bx, By, Bz, geomag_grid_parms, ...
    elvarr, azarr, maxdist, tol, refractive_ind, blind_range)

%% Raytrace_itsi
% Trace a ray in the topside, going to/from the same location
% Output:
homed_ray = NaN;

%% call raytrace to "global search" the problem
% maxdist is the maximum distance in metres
origin_lat = txloc(1);
origin_long = txloc(2);
origin_ht = txloc(3);

nhops = 1; % number of hops

rays = gs_raytrace(elvarr, azarr, freq, nhops, OX_mode, ...
    origin_lat, origin_long, origin_ht, iono_en_grid, iono_en_grid_5, ...
    collision_freq, iono_grid_parms, Bx, By, Bz, geomag_grid_parms, tol, refractive_ind);

%% Identify the closest rays
close_rays = [];
ids = [];
min_di = 1E9;

for r = 1:length(rays)
    [di, close_pt, id] = ray_dist(rays(r), rxloc);
    if di < maxdist
        rays(r).close_id = id;
        rays(r).close_pt = close_pt;
        rays(r).min_dist = di;

        close_rays = [close_rays, rays(r)];

        ids = [ids, id];
    end
    if di < min_di
        min_di = di;
    end
end

if isempty(close_rays)
    disp('No valid rays')
    return
end

% fprintf('Located %i potentially valid rays, down to %1.1f km\n', length(close_rays), min_di / 1E3)

%% Kick out those rays with range < blind_range
good_rays = [];
for r = 1:length(close_rays)
    return_range = close_rays(r).group_range(close_rays(r).close_id);
    if return_range > (blind_range / 1E3)
        good_rays = [good_rays, close_rays(r)];
    end
end

% fprintf('Located %i rays > blind range of %1.1f km\n', length(good_rays), blind_range/ 1E3)
if isempty(good_rays)
    disp('all rays within blind range')
    return
end
min_dist = nan(size(good_rays));

for r = 1:length(good_rays)
    min_dist(r) = good_rays(r).min_dist;
end
best_ray = good_rays(min_dist == min(min_dist));
best_ray = best_ray(1);

homed_ray = home_ray(best_ray, rxloc, OX_mode, nhops, tol, ...
    iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    Bx, By, Bz, geomag_grid_parms, refractive_ind);




%%
% plot_rays(txloc, rxloc, rays)

























