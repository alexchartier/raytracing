function [range, dist, path_loss_db, fspl_db] = ...
    calc_fspl(ray, loc, OX_mode, nhops, tol, ...
    iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    Bx, By, Bz, geomag_grid_parms, refractive_ind, reflect)
%% calc_fspl
% Calculate the path loss using a triad of rays surrounding the ray of
% interest, plus the free-space path loss using geometric distance and the
% standard formula. 


%% define ray triad
del = 0.1; 
elvarr = [ray.initial_elev - del, ray.initial_elev, ray.initial_elev + del];
azarr = [ray.initial_bearing - del, ray.initial_bearing, ray.initial_bearing + del];

origin_lat = loc(1); 
origin_long = loc(2); 
origin_ht = loc(3);

del_rays = gs_raytrace(elvarr, azarr, ray.frequency, nhops, OX_mode, ...
    origin_lat, origin_long, origin_ht, iono_en_grid, iono_en_grid_5, ...
    collision_freq, iono_grid_parms, Bx, By, Bz, geomag_grid_parms, tol, refractive_ind);


for r = 1:length(del_rays)
    [dl, closest_pt, id] = ray_dist(del_rays(r), loc, reflect);
    del_rays(r).close_pt = closest_pt;
    del_rays(r).min_dist = dl;
    del_rays(r).close_id = id;

    % Calculate XYZ coords of rays
    [x, y, z] = wgs84_llh2xyz(del_rays(r).lat, del_rays(r).lon, del_rays(r).height * 1E3);
    del_rays(r).x = x;
    del_rays(r).y = y;
    del_rays(r).z = z;
end

ray = del_rays(5);
del_rays = [del_rays(1), del_rays(3), del_rays(8)];


%% Poynting vector
tx_pts = nan(3, 3);
rx_pts = nan(3, 3);
small_sphere_area = nan(3, 1);
ptidx = 2;  % index of nearby ray points to consider in area calculation 

for r = 1:3
    tx_pts(1, r) = del_rays(r).x(ptidx);
    tx_pts(2, r) = del_rays(r).y(ptidx);
    tx_pts(3, r) = del_rays(r).z(ptidx);
    small_sphere_area(r) = 4 * pi * sqrt( ...
        (del_rays(r).x(ptidx) - del_rays(r).x(1)).^2 + ...
        (del_rays(r).y(ptidx) - del_rays(r).y(1)).^2 + ...
        (del_rays(r).z(ptidx) - del_rays(r).z(1)).^2) .^2;
    rx_pts(:, r) = del_rays(r).close_pt;
end

ons = [1 1 1];
Area_rx = 0.5 * sqrt(...
    det([rx_pts(1, :); rx_pts(2, :); ons])^2 + ...
    det([rx_pts(2, :); rx_pts(3, :); ons])^2 + ...
    det([rx_pts(3, :); rx_pts(1, :); ons])^2);

Area_tx= 0.5 * sqrt(...
    det([tx_pts(1, :); tx_pts(2, :); ons])^2 + ...
    det([tx_pts(2, :); tx_pts(3, :); ons])^2 + ...
    det([tx_pts(3, :); tx_pts(1, :); ons])^2);
%path_loss_db = 10 * log10(Area_tx / Area_rx);

%% Calculate the loss via free-space path loss (FSPL) and by the area increase 
S = 1 / mean(small_sphere_area) .* (Area_tx / Area_rx);
wavelength = 3E2 ./ ray.frequency;
dist = ray.geometric_distance(ray.close_id);
path_loss_db = 10 * log10(S .* wavelength.^2 / (4 .* pi));
fspl_db = -1 * fspl(dist * 1E3, wavelength);
% fprintf('free space: %1.1f dB, raytracing: %1.2f dB \n', fspl_db, path_loss_db);

%% Identify group delay
range = ray.group_range(ray.close_id);
























