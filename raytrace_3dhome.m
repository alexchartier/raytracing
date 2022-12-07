function [ray, elvarr, azarr] = raytrace_3dhome(freq, txloc, rxloc, ...
    iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    Bx, By, Bz, geomag_grid_parms, OX_mode, tol, refractive_ind, nhops, maxdist)

%% call raytrace to "global search" the problem

origin_lat = txloc(1); 
origin_long = txloc(2); 
origin_ht = txloc(3);

% elvarr = -90:3:-30;               % initial elevation of rays
% azarr = -180:3:180; % initial bearing of rays
% nhops = 1;                  % number of hops
% OX_mode = 1;
% tol(1) = 1e-7;
% tol(2) = 0.025; % km
% tol(3) = 25; % km

reflect = true; % false for transionospheric rays

% Set up the search array
az = azimuth(txloc(1), txloc(2), rxloc(1), rxloc(2));
azarr = az-3:az+3;
elvarr = 10:2:90;

rays = gs_raytrace(elvarr, azarr, freq, nhops, OX_mode, ...
    origin_lat, origin_long, origin_ht, iono_en_grid, iono_en_grid_5, ...
    collision_freq, iono_grid_parms, Bx, By, Bz, geomag_grid_parms, tol, 1);


%% Identify the closest ray
err = 1E12;
for r = 1:length(rays)
    di = ray_dist(rays(r), rxloc, reflect);
    if di < err
        err = di;
        ri = r;
    end
end

if err == 1E12
    fprintf('No valid rays\n')
    ray.home = false;
    return
end

% fprintf('Located start ray %1.1e km from receive location\n', err / 1E3)
start_ray = rays(ri);


%% Optimize from the closest ray (note start_ray contains the txloc information)

f = @(X)raytrace_err(X, start_ray, rxloc, OX_mode, nhops, tol, ...
    iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    Bx, By, Bz, geomag_grid_parms, refractive_ind, reflect);

% 'MaxFunEvals', 100, 'MaxIter', 100, 
options = optimset('TolFun', maxdist);
x0 = [start_ray.initial_elev, start_ray.initial_bearing];
[X, fval] = fminsearch(f,x0, options);

[err, ray] = f(X); % just to give back the ray


if fval < maxdist
    fprintf('Optimized to %1.1e km from receive location\n', fval/1E3)
        ray.home = true;

else
    fprintf('Failed to home\n')
    ray.home = false;
end
