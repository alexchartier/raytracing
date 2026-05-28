%% set inputs
mod_fn_fmt = '/Users/chartat1/data/sami3/2017_tid/sami_mat/{YYYY-mm-dd_HHMM}.mat';
vert_ray_fn_fmt = ['/Users/chartat1/data/sami3/2017_tid/recon_rays/O_mode/' ...
    'vert_{YYYY-mm-dd_HHMM}_%1.1fN_%1.1fE_%ikm.mat'];
oblique_ray_fn_fmt = ['/Users/chartat1/data/sami3/2017_tid/recon_rays/', ...
    'O_mode/oblique_{YYYY-mm-dd_HHMM}_%1.1fN_%1.1fE_%1.1fN_%1.1fE_%ikm.mat'];
gs_fn_fmt = ['/Users/chartat1/data/sami3/2017_tid/recon_rays_ground_space/', ...
    'O_mode/sim_twosat_rays_%s_%i.mat'];

fig1_fn_fmt = '/Users/chartat1/Documents/Papers/2024_nebula/twoball_movie/%03d.png';

modes = ['O', 'X'];
lats = 30:55;
lon = -77.5;
sc_alt = 580;
latlimit = [33, 51];
htlimit = [0, 700];

alongtrack_spacing_km = 600;
time = datenum(2017, 1, 10, 18, 0, 0);


%% Load
mod = loadstruct(filename(mod_fn_fmt, time));

dene = zeros(length(lats), length(mod.alt));

for li = 1:length(lats)
    dene(li, :) = interp_sami(mod, [lats(li), lon]);
end


%% Fig 1. two-spacecraft flyover with ground transmitter
% a: 5 km spacing
% b: 1500 km spacing
% c: 3000 km spacing

for i = 1:length(lats) - 1
    %%
    close all
    figure('units','normalized','outerposition',[0 0 1 1]);
    colormap parula
    hold on
    [~, hC] = contourf(lats, mod.alt, elec2freq(dene') / 1E3, 50);
    set(hC, 'LineStyle', 'none')
    xlabel('Lat (°)')
    ylabel('Alt (km)')
    xlim([min(lats), max(lats)])
    ylim([0, 600])
    cl = colorbar;
    ylabel(cl, 'Plasma Frequency (MHz)')
    set(gca, 'color', 'k', 'FontSize', 30)
    set(gcf, 'InvertHardCopy', 'off');
    set(gcf,'Color',[0 0 0]); % RGB values [0 0 0] indicates black color

    txloc = [lats(i), lon, sc_alt];
    refloc = [lats(i + 1), lon, sc_alt];
    rxloc = calc_rx_loc_for_ois(txloc, refloc, alongtrack_spacing_km);

    vert_rays = loadstruct(sprintf(filename(vert_ray_fn_fmt, time), txloc(1), txloc(2), txloc(3)));

    oblique_rays = loadstruct(sprintf(filename(oblique_ray_fn_fmt, time), ...
        txloc(1), txloc(2), rxloc(1), rxloc(2), txloc(3)));

    gs_rays = loadstruct(sprintf(gs_fn_fmt, 'ground_space', i));
    if i == 1
        txloc_gd = gs_rays{1}(1).txloc;
    end
    if i < length(lats) - 6
        gs_rays_2 = loadstruct(sprintf(gs_fn_fmt, 'ground_space', i+5));
        for ri = 1:length(gs_rays_2)
            for r = 1:length(gs_rays_2{ri})
                plot(gs_rays_2{ri}(r).lat, gs_rays_2{ri}(r).height, '-m', 'LineWidth', 2)
            end
        end
    end

    for ri = 1:length(vert_rays)
        for r = 1:length(vert_rays{ri})
            plot(vert_rays{ri}(r).lat, vert_rays{ri}(r).height, '-w')
        end
    end
    for ri = 1:length(oblique_rays)
        for r = 1:length(oblique_rays{ri})
            plot(oblique_rays{ri}(r).lat, oblique_rays{ri}(r).height, '-w')
        end
    end
    for ri = 1:length(gs_rays)
        for r = 1:length(gs_rays{ri})
            plot(gs_rays{ri}(r).lat, gs_rays{ri}(r).height, '-m', 'LineWidth', 2)
        end
    end

    plot(vert_rays{1}(1).txloc(1), vert_rays{1}(1).txloc(3), '.r', 'MarkerSize', 50)
    plot(oblique_rays{1}(1).rxloc(1), oblique_rays{1}(1).rxloc(3),  '.r', 'MarkerSize', 50)
    plot(txloc_gd(1), txloc_gd(3),  '.g', 'MarkerSize', 50)

    export_fig(sprintf(fig1_fn_fmt, i))
end



%% Fig 2: ionogram with O & X modes, multipath etc gets cleaned up using az/el
