%% OTHR_sim.m
% Simulate 

sami_fn_fmt = '/Users/chartat1/data/sami3/2017_tid/sami_mat/{YYYY-mm-dd_HHMM}.mat';
out_fn_fmt = '/Users/chartat1/data/raytracing/othr_depletion/rays_{YYYY-mm-dd_HHMM}.mat';
times = datenum(2017, 1, 10, 18, 0, 0):5/60/24:datenum(2017, 1, 10, 19, 10, 0);
plt_fn_fmt = '/Users/chartat1/data/raytracing/othr_depletion/plots/{YYYY-mm-dd_HHMM}.png';
image_file = 'graphics/worldmap_bw.png';

txloc = [43.23, -120.67, 0];
freq = 10;
OX_mode = 0;
nhops = 1;
tol = [1e-7 0.01 25];

elvarr = 5:40;
azarr = [270:5:360];

depletion_loc = [46, -127.5];

%% Generate rays
depletion  = [0 0 0.5 0.8 0.9 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0];
for t = 1:length(times)
    sami = loadstruct(filename(sami_fn_fmt, times(t)));

    % Add a depletion
    sami_lon = -177.5:5:177.5;
    deplete_lati = sami.lat == depletion_loc(1);
    deplete_loni = sami_lon == depletion_loc(2);

    iono_en_grid = sami.iono_en_grid;
    iono_en_grid(deplete_lati, deplete_loni, :) = ...
        iono_en_grid(deplete_lati, deplete_loni, :) * (1 - depletion(t));

    % Generate rays
    rays = gs_raytrace(elvarr, azarr, freq, nhops, OX_mode, ...
        txloc(1), txloc(2), txloc(3), iono_en_grid, iono_en_grid, ...
        sami.collision_freq, sami.iono_grid_parms, sami.Bx, sami.By, sami.Bz, ...
        sami.geomag_grid_parms, tol, 0, 0, 0, 0);


    savestruct(filename(out_fn_fmt, times(t)), rays);
    fprintf('Saved to %s\n', filename(out_fn_fmt, times(t)))    
end

%% plot
for t= 1:length(times)
    sami = loadstruct(filename(sami_fn_fmt, times(t)));

    rays = loadstruct(filename(out_fn_fmt, times(t)));
    txloc(3) = 20;
    % Ionosphere plot

    deplete_lati = sami.lat == depletion_loc(1) ;
    deplete_loni = sami_lon == depletion_loc(2);

    iono_en_grid = sami.iono_en_grid;
    iono_en_grid(deplete_lati, deplete_loni, :) = ...
        iono_en_grid(deplete_lati, deplete_loni, :) * (1 - depletion(t));

    lath = sami.lat(1):sami.lat(end);
    lonh = sami.lon(1):sami.lon(end);
    [lath3, lonh3] = meshgrid(lath, lonh);
    nmf2 = max(iono_en_grid(2:end-1, :, :), [], 3);
    nmf2_h = interp2(sami.lon, sami.lat(2:end-1), nmf2, lonh3(:), lath3(:));
    nmf2_h = reshape(nmf2_h, size(lonh3))';
    
    G = flipud(nmf2_h);

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
    view(-50, 22);

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
    zoom(4)

    camdolly(-0.5, 2, 0)

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
