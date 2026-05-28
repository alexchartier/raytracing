function [err, ray] = raytrace_err(X, start_ray, rxloc, OX_mode, nhops, tol, ...
     Ne, Bx_i, By_i, Bz_i)

%% Raytrace from Tx to Rx
% Error is difference between closest point of ray and rxloc in meters
% X = el, az
% start_ray - ray defined by pharlap as a reasonable first guess
% rxloc = lat, lon, ht

%fprintf('elv: %1.1f,  azm: %1.1f\n', X(1), X(2))

%% Extract a few params
origin_lat = start_ray.lat(1);
origin_long = start_ray.lon(1);
origin_ht = start_ray.height(1);
freq = start_ray.frequency;


%% Calculate the ray state vector
rsv = calc_rsv(origin_lat, origin_long, origin_ht, X(1), X(2), freq, ...
        OX_mode, Ne, Bx_i, By_i, Bz_i);

%% Raytrace
[ray_data, ray, ~] = ...
    raytrace_3d(origin_lat, origin_long, origin_ht, X(1), X(2), freq, ...
    OX_mode, nhops, tol, rsv);
% no need to pass in the grids as gs_raytrace has already called
% raytrace_3d by now
ray.final_Doppler = ray_data.Doppler_shift;

%%  Find minimum distance from raypath to receiver location
err = ray_dist(ray, rxloc);

% fprintf('Error: %1.1f \n', err)

