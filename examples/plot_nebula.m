%% plot the full constellation

%% set inputs
rootdir = '~/data/nebula/osse_fullsat_iri_hisolar/homings/';
ionodir = '~/data/nebula/osse_fullsat_iri_hisolar/truth/';
plotdir = '~/Documents/Papers/2024_Nebula/plots/';
R = 6371000;

%% load
filelist = dir(fullfile(rootdir, '**/*.mat'));  %get list of files and folders in any subfolder
filelist = filelist(~[filelist.isdir]);

rays = {};
for t = 1:length(filelist)
    rays{t} = loadstruct([filelist(t).folder, '/', filelist(t).name]);
end

%% Stats
n_O = 0;
n_X = 0;
dist = zeros([length(rays), 1]);
golden_fq = zeros([length(rays), 1]);
golden_dist = zeros([length(rays), 1]);
golden_ref_ind = zeros([length(rays), 1]);
golden_abs = zeros([length(rays), 1]);
golden_wnorm_B = zeros([length(rays), 1]);
golden_az = zeros([length(rays), 1]);
golden_el = zeros([length(rays), 1]);

for t = 1:length(rays)
    for r = 1:length(rays{t})
        for j = 1:length(rays{t}{r})
            % tally up the links
            if rays{t}{r}(j).OX_mode == 1
                n_O = n_O + 1;
            elseif rays{t}{r}(j).OX_mode == -1
                n_X = n_X + 1;
            end

            % get the max range
            if rays{t}{r}(j).geometric_dist_to_rx > golden_dist(t)
                sc_az = azimuth(rays{t}{r}(j).txloc(1), rays{t}{r}(j).txloc(2),...
                    rays{t}{r}(j).rxloc(1), rays{t}{r}(j).rxloc(2));
                golden_dist(t) = rays{t}{r}(j).geometric_dist_to_rx;
                golden_fq(t) = rays{t}{r}(j).frequency;
                golden_ref_ind(t) = rays{t}{r}(j).refractive_index(1);
                golden_abs(t) = rays{t}{r}(j).total_absorption;
                golden_wnorm_B(t) = min(rays{t}{r}(j).wavenorm_B_angle);
                golden_az(t) = rays{t}{r}(j).initial_bearing - sc_az ;
                golden_el(t) = rays{t}{r}(j).initial_elev;

                txloc = rays{t}{r}(j).txloc;
                rxloc = rays{t}{r}(j).rxloc;
                txloc_XYZ = sphcart([txloc(3) * 1E3 + R, deg2rad(txloc(1)), deg2rad(txloc(2))]);
                rxloc_XYZ = sphcart([rxloc(3) * 1E3 + R, deg2rad(rxloc(1)), deg2rad(rxloc(2))]);
                dist(t) = sqrt(sum((txloc_XYZ - rxloc_XYZ).^2)) / 1E3;
            end
        end
    end
end

golden_wnorm_B(golden_wnorm_B > 90) = 90 - golden_wnorm_B(golden_wnorm_B > 90);


%% histograms
close
hf = figure('Position', [10, 10, 900, 500]);

% hf = colordef(hf, 'black');
% hf.Color = 'k';

histogram2(golden_dist, golden_fq, [0:50:3200], [0:80], 'FaceColor', 'flat')
title('"Golden" frequency vs path length')
ylabel('Freq (MHz)')

y = colorbar;
ylabel(y, '# rays', 'rotation', 270, 'FontSize', 18)
view([0, 90])
grid on; grid minor
set(gca, 'FontSize', 18, 'FontName', 'Futura')
xlabel('Raypath geometric range (km)')
export_fig([plotdir, 'range_freq.png'])


%%
close
hf = figure('Position', [10, 10, 800, 1200]);

% hf = colordef(hf, 'black');
% hf.Color = 'k';

% absorption
subplot(4, 1, 1)
histogram2(dist, golden_abs, [0:50:3000], [0:30], 'FaceColor', 'flat')
title('Total Absorption')
ylabel('Absorption (dB)')

y = colorbar;
ylabel(y, '# rays', 'rotation', 270, 'FontSize', 18)
view([0, 90])
grid on; grid minor

set(gca, 'FontSize', 18, 'FontName', 'Futura', 'XtickLabels', '')

% wavenorm B
subplot(4, 1, 2)
histogram2(dist, golden_wnorm_B, [0:50:3000], [0:90], 'FaceColor', 'flat')
title('B⊥ Angle')
ylabel({'Closest angle', 'to B⊥ (°)'})
% xlabel('S/C sep (km)')

y = colorbar;
ylabel(y, '# rays', 'rotation', 270, 'FontSize', 18)
view([0, 90])
grid on; grid minor

set(gca, 'FontSize', 18, 'FontName', 'Futura', 'XTickLabel', '')

% azimuth
subplot(4, 1, 3)
histogram2(dist, golden_az, [0:50:3000], [-40:40], 'FaceColor', 'flat')
title('Received Azimuth')
ylabel({'Az. east of', 'great circle (°)'})

y = colorbar;
ylabel(y, '# rays', 'rotation', 270, 'FontSize', 18)
view([0, 90])
grid on; grid minor

set(gca, 'FontSize', 18, 'FontName', 'Futura', 'XtickLabels', '')

% Elevation
subplot(4, 1, 4)
histogram2(dist, golden_el, [0:50:3000], [-90:0], 'FaceColor', 'flat')
title('Received Elevation')
ylabel('Elevation (°)')
xlabel('S/C sep (km)')

y = colorbar;
ylabel(y, '# rays', 'rotation', 270, 'FontSize', 18)
view([0, 90])
grid on; grid minor

set(gca, 'FontSize', 18, 'FontName', 'Futura')
export_fig([plotdir, 'stats.png'])


%% Ionosphere plot
iono_en_grid = loadstruct([ionodir, 'iono_en_grid.mat']);
lat_g = loadstruct([ionodir, 'lats.mat']);
lon_g = loadstruct([ionodir, 'lons.mat']);

G = max(iono_en_grid, [], 3);
C = colormap;  % Get the figure's colormap.
L = size(C,1);
% Scale the matrix to the range of the map.
Gs = round(interp1(linspace(min(G(:)),max(G(:)),L),1:L,G));
H = reshape(C(Gs,:),[size(Gs) 3]);

space_color = 'k';
npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
alpha = 1; % globe transparency level, 1 = opaque, through 0 = invisible

% Mean spherical earth
erad    = 6371008.7714; % equatorial radius (meters)
prad    = 6371008.7714; % polar radius (meters)

figure('Color', space_color, 'Position', [10, 10, 1200, 1200]);

hold on;

% Turn off the normal axes
set(gca, 'NextPlot','add', 'Visible','off');
axis equal;
axis auto;

% Set initial view

view(0,30);

axis vis3d;

image_file = '1024px-Land_ocean_ice_2048.jpg';
[x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);
gl = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
if ~isempty(GMST0)
    hgx = hgtransform;
    set(hgx,'Matrix', makehgtform('zrotate',GMST0));
    set(gl,'Parent',hgx);
end
cdata = imread(image_file);
set(gl, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
[x, y, z] = ellipsoid(0, 0, 0, erad * 1.05, erad * 1.05, prad* 1.05, npanels);
gl2 = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 'None');
set(gl2, 'FaceColor', 'texturemap', 'CData', H, 'FaceAlpha', 0.8, 'EdgeColor', 'none');

h = colorbar('YTickLabel', linspace(round(min(G(:))), round(max(G(:))), 6));
set(h, 'Color','w', 'FontSize', 20)
export_fig([plotdir, 'iri.png'])
% close


%% Global plot with all rays
% earth_example
% for i = 1:length(rays)
%     for r = 1:length(rays{i})
%         plot_rays(rays{i}{r}, rays{i}{r}(1).txloc, rays{i}{r}(1).rxloc, 'w');
%     end
% end

%% plot the satellites

% close
dirn = '~/data/nebula/STK_positions/Nebula_Constellation_LLA_AllTimesteps/';
in_fn_fmt = [dirn, 'Nebula_Constellation_LLA_Time%i.csv'];
timestep = 10;
space_color = 'k';

% gc_lim = 2972

for t = 1% :67


    figure('Color', space_color, 'Position', [10, 10, 1200, 1200]);

    hold on;

    iono_en_grid = loadstruct([ionodir, 'iono_en_grid.mat']);
    lat_g = loadstruct([ionodir, 'lats.mat']);
    lon_g = loadstruct([ionodir, 'lons.mat']);
    % [~, hc] = contourf(lon_g, lat_g, max(iono_en_grid, [], 3), 50);
    % set(hc, 'LineStyle', 'None')
    % grid on
    % grid minor
    % xlabel('Lon (°)')
    % ylabel('Lat (°)')
    % colorbar

    G = max(iono_en_grid, [], 3);
    C = colormap;  % Get the figure's colormap.
    L = size(C,1);
    % Scale the matrix to the range of the map.
    Gs = round(interp1(linspace(min(G(:)),max(G(:)),L),1:L,G));
    H = reshape(C(Gs,:),[size(Gs) 3]);

    space_color = 'k';
    npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
    alpha = 1; % globe transparency level, 1 = opaque, through 0 = invisible

    % Mean spherical earth
    erad    = 6371008.7714; % equatorial radius (meters)
    prad    = 6371008.7714; % polar radius (meters)


    % Turn off the normal axes
    set(gca, 'NextPlot','add', 'Visible','off');
    axis equal;
    axis auto;

    % Set initial view

    view(0,30);

    axis vis3d;

    image_file = '1024px-Land_ocean_ice_2048.jpg';
    [x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);
    gl = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
    if ~isempty(GMST0)
        hgx = hgtransform;
        set(hgx,'Matrix', makehgtform('zrotate',GMST0));
        set(gl,'Parent',hgx);
    end
    cdata = imread(image_file);
    set(gl, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
    [x, y, z] = ellipsoid(0, 0, 0, erad * 1.05, erad * 1.05, prad* 1.05, npanels);
    gl2 = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 'None');
    set(gl2, 'FaceColor', 'texturemap', 'CData', H, 'FaceAlpha', 0.8, 'EdgeColor', 'none');

    h = colorbar('YTickLabel', linspace(round(min(G(:) / 1E5)) * 1E5, ...
        round(max(G(:)/1E5))*1E5, 6));
    set(h, 'Color','w', 'FontSize', 20, 'Position', [.8, .24, .03, .5])

    ylabel(h,'Peak Density (el. cm{-3})','FontSize',34,'Rotation',90);

    %%
    data = readtable(sprintf(in_fn_fmt, t));
    XYZ  = lla2ecef([data.Var2, data.Var3, data.Var4 * 1E3], 'WGS84');

    % gen rays for all < 3000km great-circle seps
    ray_ct = 0;
    for i = 1:size(XYZ, 1)
        dists = sqrt(sum((XYZ - XYZ(i, :)) .^2, 2));
        idx = dists < 3000E3;
        XYZ_i = XYZ(idx, :); %repmat(XYZ(i, :), [sum(idx), 1]);
        for j = 1:size(XYZ_i, 1)
            plot3([XYZ(i, 1), XYZ_i(j, 1)], [XYZ(i, 2), XYZ_i(j, 2)], ...
                [XYZ(i, 3), XYZ_i(j, 3)], '-w')
            ray_ct = ray_ct + 1;
        end
    end


    % plot locations

    plot3(XYZ(:, 1), XYZ(:, 2), XYZ(:, 3),'.r', 'MarkerSize', 28)


    text(0, -6E6, -6.1E6, sprintf('%i secs', (t-1) * timestep), 'color', 'w')

    % % save
    % export_fig([plotdir, sprintf('locs/%03d.png', t)])
    % close;
end


% % export_fig([plotdir, 'raycovg_5min.png'])











