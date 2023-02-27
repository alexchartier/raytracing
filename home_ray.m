function ray = home_ray(start_ray, rxloc, OX_mode, nhops, tol, ...
    iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    Bx, By, Bz, geomag_grid_parms, refractive_ind)

%% Work out if we're reflecting or not
reflect = test_reflect(start_ray.height(1), rxloc(3));

%% Optimize from the closest ray (note start_ray contains the txloc information)

f = @(X)raytrace_err(X, start_ray, rxloc, OX_mode, nhops, tol, ...
    iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    Bx, By, Bz, geomag_grid_parms, refractive_ind, reflect);

%x0 = [-90, 0];
x0 = [start_ray.initial_elev, start_ray.initial_bearing];
options = optimset('TolFun', 10);
[X, fval] = fminsearch(f, x0, options);

[err, ray] = f(X);
ray.min_dist = err;

if fval < 1E3
    fprintf('Optimized to %1.1e km from receive location\n', fval / 1E3)
        ray.home = true;
        [di, close_pt, id, group_path] = ray_dist(ray, rxloc);
        ray.group_range_to_rx = group_path;
else
    fprintf('Failed to home\n', fval)
    ray = NaN;
end