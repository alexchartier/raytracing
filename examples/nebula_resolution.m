%% nebula_resolution.m
% TODO: Calculate all midpoints generated over 10 minutes at 10-second
% resolution. Assume files are at 10-s resolution
% Get the mean and max spacing between points

% close

%% Set inputs
ionodir = '~/data/nebula/osse_fullsat_iri_hisolar/truth/';

dirn = '~/data/nebula/STK_positions/Nebula_Constellation_LLA_AllTimesteps/';
in_fn_fmt = [dirn, 'Nebula_Constellation_LLA_Time%i.csv'];

space_color = 'k';
max_spacing_m = [1, 1500, 3000] * 1E3;

latlims = 0:10:70;

Rad = 64E5;

%% Loop through maximum spacings
voronoi_dists = zeros(length(latlims) - 1, 3);
avg_density= zeros(length(latlims) - 1, 3);


for mi = 1:length(max_spacing_m)

    %% calculate midpoints
    midpts_XYZ = [];
    midpts_ll = [];
    for t = 1:60
        data = readtable(sprintf(in_fn_fmt, t));
        XYZ  = lla2ecef([data.Var2, data.Var3, data.Var4 * 1E3], 'WGS84');

        % gen rays for all < 3000km great-circle seps
        for i = 1:size(XYZ, 1)
            dists = sqrt(sum((XYZ - XYZ(i, :)) .^2, 2));
            idx = dists < max_spacing_m(mi);
            XYZ_i = XYZ(idx, :);
            SphV = cartsph(XYZ(i, :));

            lat1 = rad2deg(SphV(2));
            lon1 = rad2deg(SphV(3));

            for j = 1:size(XYZ_i, 1)
                SphV = cartsph(XYZ_i(j, :));
                lat2 = rad2deg(SphV(2));
                lon2 = rad2deg(SphV(3));
                [latMid, lonMid] = midpointLatLon(lat1, lon1, lat2, lon2);
                XYZ_mid = sphcart([SphV(1), deg2rad(latMid), deg2rad(lonMid)]);
                midpts_XYZ = [midpts_XYZ; XYZ_mid];
                midpts_ll = [midpts_ll; [latMid, lonMid]];
            end
        end
    end

    %% Calculate average spacing
    for li = 1:length(latlims) - 1
        latidx = midpts_ll(:, 1) > latlims(li) & midpts_ll(:, 1) < latlims(li + 1);

        area = areaquad(latlims(li), 0, latlims(li + 1), 360, wgs84Ellipsoid('km'));
        avg_density(li, mi) = sum(latidx) / area;
    end
    %% calculate maximum spacing
    % Create the Voronoi diagram
    [vx, vy] = voronoi(midpts_ll(:, 1), midpts_ll(:, 2));
    vx = vx(:);
    vy = vy(:);

    for li = 1:length(latlims) - 1

        latidx = midpts_ll(:, 1) > latlims(li) & midpts_ll(:, 1) < latlims(li + 1);
        idx = vx > min(midpts_ll(latidx, 1)) & vx < max(midpts_ll(latidx, 1)) & ...
            vy > min(midpts_ll(latidx, 2)) & vy < max(midpts_ll(latidx, 2));
        vx_i = vx(idx);
        vy_i = vy(idx);

        % Find the farthest Voronoi vertex
        max_distance = 0;
        center = [0 0];
        radius = 0;
        midpts_xyz = sphcart([Rad* ones(sum(latidx), 1), deg2rad(midpts_ll(latidx, :))]);

        for i = 1:length(vx_i)
            Vv_xyz = sphcart([Rad, deg2rad(vx_i(i)), deg2rad(vy_i(i))]);
            distances = sqrt(sum((Vv_xyz - midpts_xyz) .^2, 2));
            min_distance = min(distances);

            if min_distance > max_distance
                max_distance = min_distance;
                center = [vx_i(i) vy_i(i)];
                radius = max_distance;
            end
        end

        voronoi_dists(li, mi) = radius/1E3;
        disp(['Radius of the maximum inscribed circle: ' num2str(radius/1E3)]);
    end
end


%%
clf
colors = {'r', 'g', 'b'};
avg_spacing = sqrt(1 ./ (pi .* avg_density));
hold on
for i = 1:length(colors)
plot(mean([latlims(1:end-1); latlims(2:end)]), avg_spacing(:, i), ...
    ['-' colors{i}], 'LineWidth', 3)
plot(mean([latlims(1:end-1); latlims(2:end)]), voronoi_dists(:, i), ...
    ['--' colors{i}], 'LineWidth', 3)

end
hold off
xlabel('Lat (°)')
ylabel('Dist to nearest obs. (km)')
grid on
grid minor
legend({'Vertical (mean)', 'Vertical (max)', ...
    'Oblique to 1500 km (mean)', 'Oblique to 1500 km (max)', ...
    'Oblique to 3000 km (mean)', 'Oblique to 3000 km (max)'...
    })



%% plot locations
hold on
earth_example;
plot3(midpts_XYZ(:, 1), midpts_XYZ(:, 2), midpts_XYZ(:, 3),'.g', 'MarkerSize', 28)
plot3(XYZ(:, 1), XYZ(:, 2), XYZ(:, 3),'.r', 'MarkerSize', 56)






