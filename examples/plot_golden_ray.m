function plot_golden_ray(rays, title_str, fontsize)
% % Example:
% fn = '/Users/chartat1/data/nebula/iri_osse_el80-0_az0-360_2sats_hilat/homings/sat1_sat2/sim_twosat_rays_X_295.mat';
% rays = loadstruct(fn);
% title_str = fn;


arguments
    rays cell
    title_str char = ''
    fontsize double = 24
end
wgs84 = wgs84Ellipsoid("km");

colormap parula

%%  Unpack
D.ne = rays{1}.iono_en_grid;
D.lat = rays{1}.iono_lat;
D.lon = rays{1}.iono_lon;
D.alt = rays{1}.iono_alt;

rxloc = rays{1}(1).rxloc;
txloc = [rays{1}(1).lat(1), rays{1}(1).lon(1), rays{1}(1).height(1)];

%% Calculate the great-circle track
[track.lat, track.lon] = gcwaypts(txloc(1), txloc(2), rxloc(1), rxloc(2), 9);
dist = distance(txloc(1), txloc(2), rxloc(1), rxloc(2), wgs84);
track.dist = linspace(0, dist, 10);

%% interpolate
iono_en_grid = zeros(length(track.lat), length(D.alt));
for li = 1:length(track.lat)
    for hi = 1:length(D.alt)
        iono_en_grid(li, hi) = interp2(D.lon, D.lat, D.ne(:, :, hi), ...
            track.lon(li), track.lat(li), 'linear');
    end
end

%% Figure out the golden ray
maxrg = 0;
for ri = 1:length(rays)
    for r = 1:length(rays{ri})
        disp(rays{ri}(r).geometric_dist_to_rx)
        if rays{ri}(r).geometric_dist_to_rx > maxrg
            maxrg = rays{ri}(r).geometric_dist_to_rx;
            gri = ri;
            gri2 = r;
        end
    end
end

%% plot horizontal slice
hold on
contourf(track.lon, D.alt, iono_en_grid', 50, 'LineStyle', 'None');

for ri = 1:length(rays)
    for r = 1:length(rays{ri})
        plot(rays{ri}(r).lon, rays{ri}(r).height, '-w')
                if rays{ri}(r).geometric_dist_to_rx == maxrg
        plot(rays{ri}(r).lon, rays{ri}(r).height, '-y', 'LineWidth', 10)
                end

    end
end

xlabel('Lon (deg)')
ylabel('Alt (km)')
xlim([min(track.lon), max(track.lon)])

title(title_str, 'interpreter', 'none')
clim(gca, [0, max(iono_en_grid(:)) * 2]);
hC = colorbar;
hC.Label.String = 'Electron Density (el. cm^{-3})';
hC.FontSize = fontsize;
plot(txloc(2), txloc(3), '.r', 'MarkerSize', 20)
plot(rxloc(2), rxloc(3), '.g', 'MarkerSize', 20)

set(gca, 'FontSize', fontsize, 'FontName', 'Futura')
