function [err, ray] = raytrace_err(X, start_ray, rxloc, OX_mode, nhops, tol, ...
    iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    Bx, By, Bz, geomag_grid_parms, refractive_ind, reflect)

%% Raytrace from Tx to Rx
% Error is difference between closest point of ray and rxloc in meters
% X = el, az
% start_ray - ray defined by pharlap as a reasonable first guess
% rxloc = lat, lon, ht

%fprintf('elv: %1.1f,  azm: %1.1f\n', X(1), X(2))

origin_lat = start_ray.lat(1);
origin_long = start_ray.lon(1);
origin_ht = start_ray.height(1);
freq = start_ray.frequency;

if refractive_ind < 1
    [rsv.pos_x, rsv.pos_y, rsv.pos_z] = wgs84_llh2xyz(origin_lat, origin_long, origin_ht * 1E3);
    [rsv.dir_x, rsv.dir_y, rsv.dir_z] = relaz2xyz(1, X(1), X(2), origin_lat, origin_long);

    rsv.group_path = 0;
    rsv.geom_path = 0;
    rsv.phase_path = 0;
    rsv.absorption = 0;
    rsv.indep_var = 0;
    rsv.ODE_step_size = 1000;

    % scale rsv.dir_x, rsv.dir_y, rsv.dir_z by the refractive index
    rsv.dir_x =  rsv.dir_x * refractive_ind;
    rsv.dir_y =  rsv.dir_y * refractive_ind;
    rsv.dir_z =  rsv.dir_z * refractive_ind;
    [~, ray, ~] = ...
        raytrace_3d(origin_lat, origin_long, origin_ht, X(1), X(2), freq, ...
        OX_mode, nhops, tol, rsv);
else
    [~, ray, ~] = ...
        raytrace_3d(origin_lat, origin_long, origin_ht, X(1), X(2), freq, ...
        OX_mode, nhops, tol);
end

%%  Find minimum distance from raypath to receiver location
[err, closest_pt, id] = ray_dist(ray, rxloc);
%fprintf('Error: %1.1f \n', err)

