function homed_ray = home_ray_2(start_ray, rxloc, OX_mode, nhops, tol, ...
    homing_tol_m, delta_elv, delta_az, Ne, Bx_i, By_i, Bz_i)
%% home_ray_2
% - can be used in place of home_ray
%
% Ne, Bx_i, By_i, Bz_i are the local conditions used to scale rsv in case
% you are in the plasma
% iono_en_grid etc. not used as it's already in memory (you had to run the
% raytracer to get a decent first guess)
%
% Inspired by James, H. G. (2006). Effects on transionospheric HF propagation
%            observed by ISIS at middle and auroral latitudes. Advances in
%            Space Research, 38(11), 2303-2312.

%% Inputs
origin_lat = start_ray.lat(1);
origin_long = start_ray.lon(1);
origin_ht = start_ray.height(1);
freq = start_ray.frequency;
satellite_lat = rxloc(1);
satellite_lon = rxloc(2);
satellite_alt = rxloc(3);
tol(3) = 5;

az = start_ray.initial_bearing;
elev = start_ray.initial_elev;
di = start_ray.min_dist;
if di < homing_tol_m
    homed_ray = rmfield(start_ray, 'min_dist');
end

%% Optimize azimuth and elevation
ct = 0;
while di > homing_tol_m
    %fprintf('Homing err: %1.3e m  az: %1.1f el: %1.1f \n', di, az, elev)

    %% Determine the initial quad
    [elm, azm] = meshgrid(elev - delta_elv:delta_elv:elev + delta_elv, ...
        az - delta_az:delta_az:az + delta_az);
    az1d = azm(:)';
    elv1d = elm(:)';
    freqs = ones(size(az1d)) * freq;
    goodidx = [];

    for i = 1:length(elv1d)
        rsv = calc_rsv(origin_lat, origin_long, origin_ht, elv1d(i), az1d(i), freq, ...
            OX_mode, Ne, Bx_i, By_i, Bz_i);
        if isstruct(rsv)
            if isreal(rsv.dir_x + rsv.dir_y + rsv.dir_z)
                ray_state_vec_in(i) = rsv;
                goodidx = [goodidx; i];
            end
        end
    end


    %% Raytrace the quad
    [~, ray, ~] = ...
        raytrace_3d(origin_lat, origin_long, origin_ht, elv1d(goodidx), ...
        az1d(goodidx), freqs, OX_mode, nhops, tol, ray_state_vec_in(goodidx));

    dists = ones(size(ray)) * nan;
    for r = 1:length(ray)
        
        [dists(r), ~, ~, group_path, geom_path, total_absorption] = ray_dist(ray(r), rxloc);
        ray(r).group_range_to_rx = group_path;
        ray(r).geometric_dist_to_rx = geom_path;
        ray(r).total_absorption = total_absorption;
        ray(r).rxloc = rxloc;
        ray(r).perigee = min(ray(r).height);
    end


    %% Pick the best ray (if any)
    di = min(dists);
    i_min_dist = find(dists == di, 1);  % can get multiple identical values for 90° elv.
    if di < homing_tol_m
        homed_ray = ray(i_min_dist);
        continue
    end

    pplocs = pploc(ray, rxloc);
    fini = ~isnan(pplocs(:, 1));

    if sum(fini) == 0
        homed_ray = NaN;
        return
    end

    %% Consider updating az/el if any of them are better than before
    az = griddata(pplocs(fini, 1), pplocs(fini, 2), az1d(fini), rxloc(1), rxloc(2));
    elev = griddata(pplocs(fini, 1), pplocs(fini, 2), elv1d(fini), rxloc(1), rxloc(2));
    if isempty(az) || isempty(elev) || isnan(az) || isnan(elev)
        homed_ray = NaN;
        return
    end

    %% Exit in case of repeated failure
    if ct > 10
        homed_ray = NaN;
        disp('Iterated out')
        return
    end

    ct = ct + 1;
end


homed_ray.home = true;
[~, ~, ~, group_range_to_rx] = ray_dist(homed_ray, rxloc);
assert(~isnan(group_range_to_rx), 'check ray_dist interpolation')

homed_ray.group_range_to_rx = group_range_to_rx;

homed_ray.rxloc = rxloc;


%fprintf('Homed to %1.3e m\n', di)
if homed_ray.group_range_to_rx >= max(homed_ray.group_range) || ...
        isnan(homed_ray.group_range_to_rx)
        disp('Overshot somehow')
end




