function cost = iono_err_1d(X1, rays, nhops, tol, alts, lats, lons, Bx, By, Bz)
%% define cost function for retrieval
% X is a vector of length 3, used to define the ionospheric params
% Cost is a sum of the square (truth-est) distances in meters

RE = 6371E3;

%% prep 
if ~iscell(rays); rays = {rays}; end

%% Generate the ionosphere
X = repmat(X1, 1, length(lats) * length(lons));
[iono_en_grid, iono_grid_parms] = gen_ionosphere_from_coeffs(X, alts, lats, lons);

%% one-off run to set the electron density into raytrace_3d

% Local electron density
Ne = interp_3d(lats, lons, alts, iono_en_grid, rays{1}(1).txloc);

% ray state vec
rsv = calc_rsv(rays{1}(1).txloc(1), rays{1}(1).txloc(2), rays{1}(1).txloc(3), ...
    rays{1}(1).initial_elev, rays{1}(1).initial_bearing, rays{1}(1).frequency, ...
    rays{1}(1).OX_mode, Ne, ...
    rays{1}(1).geomag_x(1), rays{1}(1).geomag_y(1), rays{1}(1).geomag_z(1));

% Raytrace to put the electron density grid into memory
raytrace_3d(rays{1}(1).txloc(1), rays{1}(1).txloc(2), rays{1}(1).txloc(3), ...
    rays{1}(1).initial_elev, rays{1}(1).initial_bearing, rays{1}(1).frequency, ...
    rays{1}(1).OX_mode, nhops, tol, iono_en_grid, iono_en_grid, ...
    zeros(size(iono_en_grid)), iono_grid_parms, Bx, By, Bz, ...
    iono_grid_parms, rsv);


%% loop over 'observed' rays
dists = [];  % distance between observed and modeled ray point
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
        dists = [dists; sqrt(sum((rxloc_XYZ - loc_XYZ).^2))];

        %plot_rays_on_2d_ionosphere(rays, iono_en_grid, lats, lons, alts)

        
    end
end
fprintf('%1.1f km\n', mean(dists/1E3))


cost = sum(dists.^2);













