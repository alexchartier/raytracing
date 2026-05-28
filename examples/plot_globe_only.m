% Render a globe with opaque Earth, grey continents, white outlines, and a starfield.
% Uses sphcart from utils; add to path if needed.
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));

close all; clc;

space_color = 'k';
npanels = 180;
erad = 6371008.7714;
prad = 6371008.7714;

figure('Color', space_color, 'units','normalized','outerposition',[0 0 1 1], ...
    'Renderer', 'opengl'); % force depth buffering
ax = gca; hold(ax, 'on'); axis(ax, 'equal'); axis(ax, 'off'); axis(ax, 'vis3d');
set(ax, 'SortMethod', 'depth'); % ensure depth sorting
view(0, -80);
colormap(ax, gray(16)); % ensure neutral colormap
lighting(ax, 'none'); % disable any inherited lighting
delete(findobj(ax, 'Type', 'light'));
shading(ax, 'flat');
set(ax, 'Color', space_color); % ensure background remains opaque black
% 
% % Starfield backdrop
% nstars = 800; star_rad = erad * 10;
% dirs = randn(nstars, 3); dirs = dirs ./ vecnorm(dirs, 2, 2);
% stars = dirs * star_rad;
% scatter3(stars(:,1), stars(:,2), stars(:,3), 2, ones(nstars,3), 'filled', ...
%     'MarkerFaceAlpha', 1.0, 'MarkerEdgeAlpha', 1.0); % fully opaque stars; depth-buffered

% Opaque Earth (double shell to eliminate any bleed-through)
[x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);
surf(x, y, -z, 'FaceColor', [0 0 0], 'FaceAlpha', 1.0, 'EdgeColor', 'none', ...
    'FaceLighting', 'none', 'EdgeLighting', 'none', 'BackFaceLighting', 'unlit', ...
    'CData', [], 'CDataMode', 'manual');

hold on
[xi, yi, zi] = ellipsoid(0, 0, 0, erad * 0.995, erad * 0.995, prad * 0.995, npanels);
surf(xi, yi, -zi, 'FaceColor', [0 0 0], 'FaceAlpha', 1.0, 'EdgeColor', 'none', ...
    'FaceLighting', 'none', 'EdgeLighting', 'none', 'BackFaceLighting', 'unlit', ...
    'CData', [], 'CDataMode', 'manual');

%%
% Continents: light grey fill + white outline
try
    coast = load('coastlines');
    lat = coast.coastlat(:); lon = coast.coastlon(:);
    r = erad * 1.01; % small offset to reduce z-fighting with globe
    nanmask = isnan(lat) | isnan(lon);
    idx = 1;
    while idx <= numel(lat)
        if nanmask(idx), idx = idx + 1; continue; end
        next_nan = find(nanmask(idx:end), 1, 'first');
        if isempty(next_nan)
            seg_end = numel(lat);
        else
            seg_end = idx + next_nan - 2;
        end
        lat_seg = lat(idx:seg_end); lon_seg = lon(idx:seg_end);
        if numel(lat_seg) >= 3
            lat_seg(end+1) = lat_seg(1);
            lon_seg(end+1) = lon_seg(1);
            seg_xyz = sphcart([r * ones(numel(lat_seg), 1), deg2rad(lat_seg), deg2rad(lon_seg)]);
            patch('XData', seg_xyz(:,1), 'YData', seg_xyz(:,2), 'ZData', -seg_xyz(:,3), ...
                'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.35, 'EdgeColor', 'none', ...
                'FaceLighting', 'none', 'EdgeLighting', 'none'); % keep continent fill flat to avoid shading artifacts
            plot3(seg_xyz(:,1), seg_xyz(:,2), -seg_xyz(:,3), 'w', 'LineWidth', 1.0);
        end
        idx = seg_end + 1;
    end
catch
    % fallback: equator outline
    th = linspace(0, 2*pi, 360);
    eq = [erad*cos(th'); erad*sin(th'); zeros(size(th'))]';
    plot3(eq(:,1), eq(:,2), -eq(:,3), 'w', 'LineWidth', 1.2);
end

shading flat
