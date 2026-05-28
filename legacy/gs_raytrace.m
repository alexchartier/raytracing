function rays = gs_raytrace(elvarr, azarr, freq, nhops, OX_mode, ...
    origin_lat, origin_long, origin_ht, iono_en_grid, iono_en_grid_5, ...
    collision_freq, iono_grid_parms, Bx, By, Bz, geomag_grid_parms, tol, ...
    refractive_ind)

%% Global search raytracing 

[elv2d, az2d] = meshgrid(elvarr, azarr);
elv1d = double(elv2d(:))';
az1d = double(az2d(:))';
freqs = ones(size(elv1d))*freq;

%% Create starting raypath state vector
clear ray_state_vec_in
for i = 1:length(elv1d)
    [rsv.pos_x, rsv.pos_y, rsv.pos_z] = wgs84_llh2xyz(origin_lat, origin_long, origin_ht * 1E3);
    [rsv.dir_x, rsv.dir_y, rsv.dir_z] = relaz2xyz(1, elv1d(i), az1d(i), origin_lat, origin_long);
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

    ray_state_vec_in(i) = rsv;
end

% Generate rays
[ray_data, rays, ray_state_vec] = ...
    raytrace_3d(origin_lat, origin_long, origin_ht, elv1d, az1d, freqs, ...
        OX_mode, nhops, tol, iono_en_grid, iono_en_grid_5, ...
        collision_freq, iono_grid_parms, Bx, By, Bz, ...
        geomag_grid_parms, ray_state_vec_in);


    
