%% sporadic E synth.m
% Simulate two satellite topside sounder mission data with a groundbased
% transmitter there too
% see also: gen_iono_grid_sporadice.m

clear

%% inputs
time = datenum(2012, 1, 1, 0, 0, 0);
mod_fn_fmt = '/Users/chartat1/data/nebula/model_input/iri_es/{YYYY-mm-dd_HHMM}.mat';
oblique_ray_fn_fmt = '/Users/chartat1/data/nebula/rays/sporadice_sim/oblique_rays_%i.mat';
vert_ray_fn_fmt = '/Users/chartat1/data/nebula/rays/sporadice_sim/vert_rays_%i.mat';
fig1_fn_fmt = '/Users/chartat1/Documents/Papers/2024_nebula/sporadice_movie/%03d.png';


freqs = 2:0.1:20;

sat_alt = 580;
satlats = 30:50;
satlon = 15;

sat_spacing = 5;

OX_mode = 1;

maxdist = 1E8;  % meters from homing
blind_range = 30E3; % two-way (e.g. there-and-back range)
tol = [1e-7 0.01 15];       % ODE solver tolerance and min max stepsizes
homing_tol_m = 1500;

B_ht_inc = 10;
B_lat_inc = 2;
B_lon_inc = 2;


%% sat-to-sat raytracing example
clear oblique_rays
fprintf('Started a new raytracing expt\n')
model = loadstruct(filename(mod_fn_fmt, time));
for l = 1:length(satlats) - sat_spacing
    txloc = [satlats(l), satlon, sat_alt];
    rxloc = [satlats(l) + sat_spacing, satlon, sat_alt];
    oblique_rays = gen_ionogram(freqs, OX_mode, txloc, rxloc, ...
        model.iono_en_grid, model.iono_en_grid, model.collision_freq, ...
        model.iono_grid_parms, model.geomag_grid_parms, ...
        model.Bx, model.By, model.Bz, ...
        maxdist, tol, homing_tol_m, true, false);
    savestruct(sprintf(oblique_ray_fn_fmt, l), oblique_rays)

    vert_rays = gen_ionogram(freqs, OX_mode, txloc, txloc, ...
        model.iono_en_grid, model.iono_en_grid, model.collision_freq, ...
        model.iono_grid_parms, model.geomag_grid_parms, ...
        model.Bx, model.By, model.Bz, ...
        maxdist, tol, homing_tol_m);
    savestruct(sprintf(vert_ray_fn_fmt, l), vert_rays)
end
% plot_rays(sat_to_sat_rays, txloc, rxloc)


%% Plotting ionospheric density slice

dene = model.dene(:, :, 1)';
for l = 1:length(satlats) - 5
    close all
    figure('units','normalized','outerposition',[0 0 1 1]);
    hold on
    subplot(2, 2, [2, 4])
    colormap parula
    hold on
    [~, hC] = contourf(model.lat, model.alt, elec2freq(dene') / 1E3, 50);
    set(hC, 'LineStyle', 'none')
    xlabel('Lat (°)')
    ylabel('Alt (km)')
    xlim([min(satlats), max(satlats)])
    ylim([0, 600])
    cl = colorbar;
    ylabel(cl, 'Plasma Frequency (MHz)')

    set(gca, 'color', 'k', 'FontSize', 30, 'XColor', 'w', 'YColor', 'w')
    set(cl, 'XColor', 'w', 'YColor', 'w')

    set(gcf, 'InvertHardCopy', 'off');
    set(gcf,'Color',[0 0 0]); % RGB values [0 0 0] indicates black color


    vert_rays = loadstruct(sprintf(vert_ray_fn_fmt, l));
    oblique_rays = loadstruct(sprintf(oblique_ray_fn_fmt, l));

    for ri = 1:length(vert_rays)
        for r = 1:length(vert_rays{ri})
            if vert_rays{ri}(r).group_range_to_rx > 2000
                continue
            end
            plot(vert_rays{ri}(r).lat, vert_rays{ri}(r).height, '-w')
        end
    end
    for ri = 1:length(oblique_rays)
        for r = 1:length(oblique_rays{ri})
            plot(oblique_rays{ri}(r).lat, oblique_rays{ri}(r).height, '-w')
        end
    end

    plot(vert_rays{1}(1).txloc(1), vert_rays{1}(1).txloc(3), '.r', 'MarkerSize', 50)
    plot(oblique_rays{1}(1).rxloc(1), oblique_rays{1}(1).rxloc(3),  '.r', 'MarkerSize', 50)

    % Vertical
    subplot(2, 2, 1)
    hold on
    for ri = 1:length(vert_rays)
        for r = 1:length(vert_rays{ri})
            plot(vert_rays{ri}(r).frequency, vert_rays{ri}(r).group_range_to_rx, ...
                '.w', 'MarkerSize', 10)
        end
    end 
    title('Vertical', 'color', 'w')
    set(gca, 'color', 'k', 'FontSize', 30, 'XColor', 'w', 'YColor', 'w', ...
        'Xlim', [2, 10], 'YLim', [0, 2000])
    xlabel('Freq (MHz)')
    ylabel('Virtual Range (km)')
    set(gca, 'YDir','reverse')
    grid on
    grid minor

    % Oblique
    subplot(2, 2, 3)
    hold on
    for ri = 1:length(oblique_rays)
        for r = 1:length(oblique_rays{ri})
            plot(oblique_rays{ri}(r).frequency, oblique_rays{ri}(r).group_range_to_rx, ...
                '.w', 'MarkerSize', 10)
        end
    end
    title('Oblique', 'color', 'w')
    set(gca, 'color', 'k', 'FontSize', 30, 'XColor', 'w', 'YColor', 'w', ...
        'Xlim', [2, 12], 'YLim', [0, 2000])
    xlabel('Freq (MHz)')
    ylabel('Virtual Range (km)')
    set(gca, 'YDir','reverse')
    grid on
    grid minor

    export_fig(sprintf(fig1_fn_fmt, l))

end
