%% OTHR_sim.m
% Simulate 

sami_fn_fmt = '/Users/chartat1/data/sami3/2017_tid/sami_mat/{YYYY-mm-dd_HHMM}.mat';
out_fn_fmt = '/Users/chartat1/data/raytracing/othr/rays/rays_{YYYY-mm-dd_HHMM}.mat';
times = datenum(2017, 1, 10, 12, 5, 0):15/60/24:datenum(2017, 1, 10, 13, 0, 0);
plt_fn_fmt = '/Users/chartat1/data/raytracing/othr/plots/{YYYY-mm-dd_HHMM}.png';


txloc = [43.23, -120.67, 0];
freq = 10;
OX_mode = 0;
nhops = 1;
tol = [1e-7 0.01 25];

elvarr = 5:40;
azarr = [270:5:360];

%% Generate rays
for t = 1:length(times)
    sami = loadstruct(filename(sami_fn_fmt, times(t)));
    
    % Generate rays
    rays = gs_raytrace(elvarr, azarr, freq, nhops, OX_mode, ...
        txloc(1), txloc(2), txloc(3), sami.iono_en_grid, sami.iono_en_grid, ...
        sami.collision_freq, sami.iono_grid_parms, sami.Bx, sami.By, sami.Bz, ...
        sami.geomag_grid_parms, tol, 0, 0, 0, 0);

    savestruct(filename(out_fn_fmt, times(t)), rays);
    fprintf('Saved to %s\n', filename(out_fn_fmt, times(t)))
end

%% plot
for t= 1% :length(times)
    sami = loadstruct(filename(sami_fn_fmt, times(t)));

    rays = loadstruct(filename(out_fn_fmt, times(t)));
    txloc(3) = 20;
    % Ionosphere plot
    
    G = flipud(max(sami.iono_en_grid(2:end-1, :, :), [], 3));

    % G = squeeze(max(sami.dene, [], 1));
    G(isnan(G)) = 1;
    G = elec2freq(G) * 1E3;
    colormap("jet")
    C = colormap;  % Get the figure's colormap.
    L = size(C,1);

    % Scale the matrix to the range of the map.
    Gs = round(interp1(linspace(min(G(:)), max(G(:)), L), 1:L, G));
    Gs(isnan(Gs)) = 1;
    H = reshape(C(Gs(:), :), [size(Gs) 3]);


    space_color = 'k';
    npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
    alpha = 1; % globe transparency level, 1 = opaque, through 0 = invisible

    % Mean spherical earth
    erad    = 6371008.7714; % equatorial radius (meters)
    prad    = 6371008.7714; % polar radius (meters)
    close

    figure('Color', space_color, 'units','normalized','outerposition',[0 0 1 1]);

    hold on;

    % Turn off the normal axes
    set(gca, 'NextPlot','add', 'Visible','off');
    axis equal;
    axis auto;

    % Set initial view

    axis vis3d;
    view(-50, 15);

    % image_file = '1024px-Land_ocean_ice_2048.jpg';
    image_file = 'worldmap_bw.png';
    [x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);
    gl = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0*[1 1 1]);

    colormap gray
    cdata = imread(image_file);
    set(gl, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
    [x, y, z] = ellipsoid(0, 0, 0, erad * 1.05, erad * 1.05, prad* 1.05, npanels);

    % colormap parula
    gl2 = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 'None');
    set(gl2, 'FaceColor', 'texturemap', 'CData', H, 'FaceAlpha', 0.65, 'EdgeColor', 'none');
    ax = gca;
    ax.Clipping = "off";
    zoom on
    zoom(3)

    camdolly(-0.5, 2, 0)

    % colormap("parula")
    % h = colorbar('YTickLabel', linspace(round(min(G(:)/1E6)), round(max(G(:)/1E6)), 6));
    %
    % ylabel(h, 'FoF2 (MHz)')
    % set(h, 'Color','w', 'FontSize', 20)

    plot_rays(rays, txloc, txloc, 'm')
    hold on
    Re = 6380E3;
    for r = 1:length(rays)

        if rays(r).height(end) <=10

            cart = sphcart([rays(r).height(end) * 1E3 + Re, ...
                deg2rad(rays(r).lat(end)), ...
                deg2rad(rays(r).lon(end))]);
            h3 = plot3(cart(1), cart(2), cart(3), 'ro', 'markersize', 10, 'markerfacecolor', 'm');
        end
    end
    hold off


    text(-8000E3, 0, Re + 2000E3, filename('{yyyy-mm-dd HH:MM} UT', times(t)), ...
        'color', 'w', 'FontSize', 40)

    export_fig(filename(plt_fn_fmt, times(t)))
    pause(0.1)
    %%
    % clf

end



%% sanity check the grid
nmax = zeros(size(times));
nmax_mod = zeros(size(times));

nmax_mod_orig = zeros(size(times));

for t = 1:length(times)
    sami = loadstruct(filename(sami_fn_fmt, times(t)));
        rays = loadstruct(filename(out_fn_fmt, times(t)));
    lat = arange(sami.iono_grid_parms(1), sami.iono_grid_parms(2), sami.iono_grid_parms(3));
    lon = arange(sami.iono_grid_parms(4), sami.iono_grid_parms(5), sami.iono_grid_parms(6));
    G = max(sami.iono_en_grid(2:end-1, :, :), [], 3);
    G(isnan(G)) = 1;
    nmax_mod(t) = interp2(lon, lat(2:end-1), G, txloc(2), txloc(1));
    nmax_mod_orig(t) = interp2(sami.lon, sami.lat, squeeze(max(sami.dene, [], 1)), txloc(2) + 360, txloc(1));
    
    for r = 1:length(rays)
        if max(rays(r).electron_density) > nmax(t)
            nmax(t) = max(rays(r).electron_density);
        end
    end

end


hold on
plot(times, nmax_mod, 'k')
plot(times, nmax, 'r')
plot(times, nmax_mod_orig, 'm')

datetick
grid on
grid minor
hold off


%% 

t = 91;

    sami = loadstruct(filename(sami_fn_fmt, times(t)));
   lat = arange(sami.iono_grid_parms(1), sami.iono_grid_parms(2), sami.iono_grid_parms(3));
    lon = arange(sami.iono_grid_parms(4), sami.iono_grid_parms(5), sami.iono_grid_parms(6));
        G = max(sami.iono_en_grid(2:end-1, :, :), [], 3);
    G(isnan(G)) = 1;
hold on
    [~, hC] = contourf(lon, lat(2:end-1), elec2freq(G) / 1E3, 50);
    
    set(hC, 'LineStyle', 'none')
    plot(txloc(2), txloc(1), '.m', 'MarkerSize', 30)
    colorbar
    colormap jet
    title(filename('{yyyy-mm-dd HH:MM} UT', times(t)))




