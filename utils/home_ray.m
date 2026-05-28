function ray = home_ray(start_ray, rxloc, OX_mode, nhops, tol, homing_tol_m, ...
     Ne, Bx_i, By_i, Bz_i, varargin)
%% Home a ray to the rxloc

% %% Work out if we're reflecting or not
% reflect = test_reflect(start_ray.height(1), rxloc(3));

% homing_tol_m = 500;


%% Set optimization options
if length(varargin) > 0
    options = varargin{1};
else
    options = optimset('TolFun', homing_tol_m, 'MaxFunEvals', 100, ...
        'TolX', 0.005, 'Display', 'off');
end


%% Optimize from the closest ray (note start_ray contains the txloc information)

f = @(X)raytrace_err(X, start_ray, rxloc, OX_mode, nhops, tol, Ne, Bx_i, By_i, Bz_i);
x0 = [start_ray.initial_elev, start_ray.initial_bearing];

% safer to remove maxfunevals and tolx, but could set to 200/0.001 for same
% effect
[X, fval] = fminsearch(f, x0, options);

[err, ray] = f(X);
ray.min_dist = err;

if err < homing_tol_m  % NOTE: was 1E3 before
%     fprintf('Optimized to %1.1e km from receive location\n', fval / 1E3)
        ray.home = true;
        [~, ~, ~, group_path, geom_path, total_absorption] = ray_dist(ray, rxloc);
        ray.group_range_to_rx = group_path;
        ray.geometric_dist_to_rx = geom_path;
        ray.total_absorption = total_absorption;
        ray.rxloc = rxloc;
        ray.perigee = min(ray.height);
else
%     fprintf('Failed to home\n', fval)
    ray = NaN;
end