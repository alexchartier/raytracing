function plot_rays_on_2d_ionosphere(rays, iono_en_grid, lats, lons, alts, title_str, fontsize)
% % Example:
% fn = '/Users/chartat1/data/nebula/iri_osse_el80-0_az0-360_2sats_hilat/homings/sat1_sat2/sim_twosat_rays_X_295.mat';
% rays = loadstruct(fn);
% title_str = fn;
% plot_rays_on_2d_ionosphere(rays, iono_en_grid, lats, lons, alts)
arguments
    rays cell
    iono_en_grid (:, :, :) double
    lats double
    lons double
    alts double
    title_str char = ''
    fontsize double = 24
end
wgs84 = wgs84Ellipsoid("km");

colormap parula

%%  Unpack
% D.ne = rays{1}.iono_en_grid;
% D.lat = rays{1}.iono_lat;
% D.lon = rays{1}.iono_lon;
% D.alt = rays{1}.iono_alt;

rxloc = rays{1}(1).rxloc;
txloc = [rays{1}(1).lat(1), rays{1}(1).lon(1), rays{1}(1).height(1)];



%% Calculate the great-circle track
[track.lat, track.lon] = gcwaypts(txloc(1), txloc(2), rxloc(1), rxloc(2), 9);
dist = distance(txloc(1), txloc(2), rxloc(1), rxloc(2), wgs84);
track.dist = linspace(0, dist, 10);

%% interpolate electron density to track
Ne = zeros(length(track.lat), length(alts));
for li = 1:length(track.lat)
    for hi = 1:length(alts)
        Ne(li, hi) = interp2(lons, lats, iono_en_grid(:, :, hi), ...
            track.lon(li), track.lat(li), 'linear');
    end
end

%% interpolate rays to track
int_rays = {};
for ri = 1:length(rays)
    for r = 1:length(rays{ri})

        if length(unique(track.lon)) == 1
            rays{ri}(r).track_alt = interp1(...
                rays{ri}(r).lat, rays{ri}(r).height, track.lat);
        elseif length(unique(track.lat)) == 1
            rays{ri}(r).track_alt = interp1(...
                rays{ri}(r).lon, rays{ri}(r).height, track.lon);
        else
            rays{ri}(r).track_alt = interp1(...
                rays{ri}(r).lat, rays{ri}(r).height, ...
                track.lat);
        end
    end
end

%% plot horizontal slice
if all(round(txloc(1:2) * 1E3) == round(rxloc(1:2) * 1E3))  % vertical mode, just do lat/alt slice
    subplot(2, 1, 1)
    iono_slice = zeros([length(lats), length(alts)]);
    for i =1:length(lats)
        for j = 1:length(alts)
            iono_slice(i, j) = interp1(lons, iono_en_grid(i, :, j), rxloc(2));
        end
    end

    hold on
    contourf(lats, alts, iono_slice', 50, 'LineStyle', 'None');

    lmin = 1E3;
    lmax = 0;
    for ri = 1:length(rays)
        for r = 1:length(rays{ri})
            hi = rays{ri}(r).height <= txloc(3);
            plot(rays{ri}(r).lat(hi), rays{ri}(r).height(hi), '-w', 'LineWidth', 3)
            if lmin > min(rays{ri}(r).lat)
                lmin = min(rays{ri}(r).lat);
            end
            if lmax < max(rays{ri}(r).lat)
                lmax = max(rays{ri}(r).lat);
            end
        end
    end

    plot(txloc(1), txloc(3), '.r', 'MarkerSize', 20)
    xlim([lmin - 3, lmax + 3])
    ylim([0, txloc(3)])
    xlabel('Lat (deg)')
    ylabel('Alt (km)')
    
    hC = colorbar;
    hC.Label.String = 'Density (el. cm^{-3})';
    hC.FontSize = fontsize;
    set(gca, 'FontSize', fontsize, 'FontName', 'Futura')
    clim(gca, [0, max(Ne(:)) * 1.1]);
   

else % oblique mode, plot along-track vs alt
    subplot(2, 1, 1)
    hold on
    contourf(track.dist, alts, Ne', 50, 'LineStyle', 'None');

    for ri = 1:length(rays)
        for r = 1:length(rays{ri})
            plot(track.dist, rays{ri}(r).track_alt, '-w', 'LineWidth', 3)
        end
    end

    xlabel('Along-track dist (km)')
    ylabel('Alt (km)')
    % xlim([min(track.lon), max(track.lon)])

    title(title_str, 'interpreter', 'none')
    clim(gca, [0, max(Ne(:))]);
    hC = colorbar;
    hC.Label.String = 'Electron Density (el. cm^{-3})';
    hC.FontSize = fontsize;
    plot(0, txloc(3), '.r', 'MarkerSize', 20)
    plot(max(track.dist), rxloc(3), '.g', 'MarkerSize', 20)

    set(gca, 'FontSize', fontsize, 'FontName', 'Futura')
end

%% plot from above (lat/lon over nmf2)
subplot(2, 1, 2)
nmf2 = max(iono_en_grid, [], 3);
lon_step = unique(diff(lons));
lat_step = unique(diff(lats));
lonlim = [min(track.lon) - lon_step, max(track.lon) + lon_step];
latlim = [min(track.lat) - lat_step, max(track.lat) + lat_step];

loni = lons >= floor(lonlim(1)) & lons <= ceil(lonlim(2));
lati = lats >= floor(latlim(1)) & lats <= ceil(latlim(2));

hold on
contourf(lons(loni), lats(lati), nmf2(lati, loni), 50, 'LineStyle', 'None');

for ri = 1:length(rays)
    for r = 1:length(rays{ri})
        hi = rays{ri}(r).height <= txloc(3);
        plot(rays{ri}(r).lon(hi), rays{ri}(r).lat(hi), '-w', 'LineWidth', 3)
    end
end

plot(txloc(2), txloc(1), '.r', 'MarkerSize', 20)
plot(rxloc(2), rxloc(1), '.g', 'MarkerSize', 20)
xlabel('Lon (deg)')
ylabel('Lat (deg)')
%xlim(lonlim)
%ylim(latlim)
% clim(gca, [0, max(nmf2(:))]);
hC = colorbar;
hC.Label.String = 'Peak Density (el. cm^{-3})';
hC.FontSize = fontsize;
set(gca, 'FontSize', fontsize, 'FontName', 'Futura')



