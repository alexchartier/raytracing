function rays = gs_raytrace(elvarr, azarr, freq, nhops, OX_mode, ...
    origin_lat, origin_long, origin_ht, iono_en_grid, iono_en_grid_5, ...
    collision_freq, iono_grid_parms, Bx, By, Bz, geomag_grid_parms, tol, ...
    Ne, Bx_i, By_i, Bz_i)

%% Global search raytracing 
% [lats, lons, alts] = parms_to_lla(iono_grid_parms);  % lla of Ne
% 
% if txloc(3) < min(alts)
%     Ne = 0;
%     Bx_i = 0;
%     By_i = 0;
%     Bz_i = 0;
% else
%     Ne = interp_3d(lats, lons, alts, iono_en_grid, txloc);
%     assert(~isnan(Ne), 'Interpolation failed')
%     [B_lat, B_lon, B_ht] = parms_to_lla(geomag_grid_parms);  % lla of B
%     Bx_i = interp_3d(B_lat, B_lon, B_ht, Bx, txloc);
%     By_i = interp_3d(B_lat, B_lon, B_ht, By, txloc);
%     Bz_i = interp_3d(B_lat, B_lon, B_ht, Bz, txloc);
% end
%
% rays = gs_raytrace(elvarr, azarr, freq, nhops, OX_mode, ...
%     origin_lat, origin_long, origin_ht, iono_en_grid, iono_en_grid_5, ...
%     collision_freq, iono_grid_parms, Bx, By, Bz, geomag_grid_parms, tol, ...
%     Ne, Bx_i, By_i, Bz_i)

%%
[elv2d, az2d] = meshgrid(elvarr, azarr);
elv1d = double(elv2d(:))';
az1d = double(az2d(:))';
freqs = ones(size(elv1d)) .* freq;


%% Create starting raypath state vector
clear ray_state_vec_in
goodidx = [];
for i = 1:length(elv1d)
    rsv = calc_rsv(origin_lat, origin_long, origin_ht, elv1d(i), az1d(i), ...
        freq, OX_mode, Ne, Bx_i, By_i, Bz_i);
    if isstruct(rsv)
        if ~isnan(rsv.dir_x + rsv.dir_y + rsv.dir_z) &&...
                isreal(rsv.dir_x + rsv.dir_y + rsv.dir_z)
            ray_state_vec_in(i) = rsv;
            goodidx = [goodidx; i];
        end
    end
end

if ~exist('ray_state_vec_in', 'var')  % rsv only returns a struct if refractive_ind > 0
    rays = NaN;
    return
end

ray_state_vec_in = ray_state_vec_in(goodidx);
elv1d = elv1d(goodidx);
az1d = az1d(goodidx);
freqs = freqs(goodidx);

% Generate rays
[~, rays, ~] = ...
    raytrace_3d(origin_lat, origin_long, origin_ht, elv1d, az1d, freqs, ...
        OX_mode, nhops, tol, iono_en_grid, iono_en_grid_5, ...
        collision_freq, iono_grid_parms, Bx, By, Bz, ...
        geomag_grid_parms, ray_state_vec_in);








