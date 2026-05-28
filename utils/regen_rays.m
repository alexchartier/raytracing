function [output_rays, dists] = regen_rays(rays, alts, lats, lons, iono_en_grid, nhops, tol)
%% (re)generate rays from an ionosphere and previous ray array 
% (potentially from a different ionosphere)

RE = 6371E3;
homing_tol_m = 3000;
%%
dists = [];  % distance between observed and modeled ray point
output_rays = cell(length(rays), 1);
for ri = 1:length(rays)
    for r = 1:length(rays{ri})
        % Local electron density
        Ne = interp_3d(lats, lons, alts, iono_en_grid, rays{ri}(r).txloc);
        
        %% Calculate the ray state vector
        rsv = calc_rsv(rays{ri}(r).txloc(1), rays{ri}(r).txloc(2), rays{ri}(r).txloc(3), ...
            rays{ri}(r).initial_elev, rays{ri}(r).initial_bearing, rays{ri}(r).frequency, ...
            rays{ri}(r).OX_mode, Ne, ...
            rays{ri}(r).geomag_x(1), rays{ri}(r).geomag_y(1), rays{ri}(r).geomag_z(1));


        %% Raytrace through model ionosphere
        [~, ray, ~] = ...
            raytrace_3d(rays{1}(1).txloc(1), rays{1}(1).txloc(2), rays{1}(1).txloc(3), ...
            rays{1}(1).initial_elev, rays{1}(1).initial_bearing, rays{1}(1).frequency, ...
            rays{1}(1).OX_mode, nhops, tol, rsv);
        ray.OX_mode = rays{1}(1).OX_mode;

        %% go through and get the modeled ray point that match the observed range
        % group range corresponds to the time-of-flight measurement
        fini = ~isnan(ray.lat);
        ray_XYZ = sphcart([ray.height(fini) * 1E3 + RE; ...
            deg2rad(ray.lat(fini)); deg2rad(ray.lon(fini))]');

        loc_XYZ = zeros(1, 3); % Location of the point along the ray
        for i = 1:3
            loc_XYZ(i) = interp1(ray.group_range(fini), ray_XYZ(:, i), ...
                rays{ri}(r).group_range_to_rx, 'linear', 'extrap');
        end

        rxloc = rays{ri}(r).rxloc;
        % %% plotting
        % loc_sph = cartsph(loc_XYZ');
        % loc_lla = [rad2deg(loc_sph(2)), rad2deg(loc_sph(3)), loc_sph(1)/1E3 - 6371];
        % earth_example
        % plot_rays(ray, loc_lla, rays{ri}(r).rxloc)

        %% Compare against observed location in XYZ
        rxloc_XYZ = sphcart([rxloc(3) * 1E3 + RE, deg2rad(rxloc(1)), deg2rad(rxloc(2))]);
        dist = sqrt(sum((rxloc_XYZ - loc_XYZ).^2));
        dists = [dists; dist];

        %plot_rays_on_2d_ionosphere(rays, iono_en_grid, lats, lons, alts)



        %% Calculate homing distance etc
        % if dist < homing_tol_m  % NOTE: was 1E3 before
        %     %     fprintf('Optimized to %1.1e km from receive location\n', fval / 1E3)
        %     ray.home = true;
        %     
        % end
        ray.group_range_to_rx = rays{ri}(r).group_range_to_rx;
        ray.perigee = min(ray.height);
        try
        output_rays{ri}(r) = ray;
        catch
            disp(1)
        end
        
    end
end