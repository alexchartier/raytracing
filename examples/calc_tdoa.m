%% calc_tdoa.m
% Calculate the time-difference of arrival through the 'truth' model and
% through various reconstructed ionospheres

% TODO: 
% - save out figures
% - See why IRI is not looking worse (looks like rays should not pass
% straight through it, maybe iono_en_grid is not matching mod.dene somehow?)

%% Set inputs
mod_fn_fmt = '/Users/chartat1/data/sami3/2017_tid/sami_mat/{YYYY-mm-dd_HHMM}.mat';
vert_recon_fn_fmt = ['/Users/chartat1/data/sami3/2017_tid/recon_iono/', ...
    'vert_recon_{YYYY-mm-dd_HHMM}_%1.1fN_%1.1fE.mat'];
oblique_recon_fn_fmt = ['/Users/chartat1/data/sami3/2017_tid/recon_iono/', ...
    'oblique_recon_{YYYY-mm-dd_HHMM}_%1.0fN_%1.0fE.mat'];
iri_fn_fmt = '/Users/chartat1/data/iri/{YYYY-mm-dd_HHMM}.mat';

time = datenum(2017, 1, 10, 18, 0, 0);
reflat = 44;  % reference latitude for the oblique reconstruction
reflat_vert = 40.5; % reference latitude for the vertical reconstruction
reflon = -77; % reference longitude for the reconstruction

txloc = [38, -77.5, 0];
rxloc_1 = [35, -77.5, 400];
rxloc_2 = [45, -77.5, 400];
freq = 10;
RE = 6371E3;

%% Set up inputs to the TDOA calculation
modlist = {...
    'Truth', ...
    'T+5', ...
    'Oblique', ...
    'Vertical', ...
    'IRI' ...
    };

flist = {...
    filename(mod_fn_fmt, time), ...
    filename(mod_fn_fmt, time + 5/60/24), ...
    filename(sprintf(oblique_recon_fn_fmt, reflat, reflon), time), ...
    filename(sprintf(vert_recon_fn_fmt, reflat_vert, reflon), time), ...
    filename(iri_fn_fmt, time) ...
    };

txloc_xyz = sphcart([RE + txloc(3), deg2rad(txloc(1)), deg2rad(txloc(2))]);
rxloc_xyz = [sphcart([RE + rxloc_1(3) * 1E3, deg2rad(rxloc_1(1)), deg2rad(rxloc_1(2))]); ...
    sphcart([RE + rxloc_2(3) * 1E3, deg2rad(rxloc_2(1)), deg2rad(rxloc_2(2))])];


%% Generate delays
delays = zeros(size(rxloc_xyz, 1), length(flist));
for i = 1:length(flist)
    fprintf('\n\n###### %s ######\n', modlist{i})

    for ri = 1:size(rxloc_xyz, 1)
        [delay_mat, rays] = calc_ionospheric_delay(flist{i}, txloc_xyz, rxloc_xyz(ri, :), freq);
        delays(ri, i) = min(delay_mat(:));
        homed_rays{i, ri} = rays;

    end
end

%% report out
TDOAs = diff(delays, [], 1);
TDOAs = TDOAs - TDOAs(1);
for i = 1:length(flist)
    fprintf('%s: %1.1f us error\n', modlist{i}, TDOAs(i))
end

%% plot the different reconstructions

%% Contour plotting
% clf
% subplot(1, 2, 1)
% contourf(sami.lat, sami.alt, elec2freq(sami.dene(:, :, sami.lon == 282.5)) /1E3, 50)
% title('Model Truth')
% clim([0 6])
% xlim([30, 55])
% xlabel('Lat (°)')
% colorbar

lats = 33:51;
lon = -77.5;
latlimit = [33, 51];
htlimit = [0, 400];
close
figure('units','normalized','outerposition',[0 0 1 1]);
hmf2 = zeros(length(flist), length(lats));
nmf2 = zeros(length(flist), length(lats));

for i = 1:length(flist)

    mod = loadstruct(flist{i});
    if i == 5
        mod.dene = cat(3, mod.dene(:, :, mod.lon > 180), mod.dene(:, :, mod.lon <= 180));
    end
    mod.lon(mod.lon < 0) = mod.lon(mod.lon < 0) + 360;

    % get the density onto the required track
    % mod.lon(mod.lon > 180) = mod.lon(mod.lon > 180) - 360;
    % loni = mod.lon == txloc(2);
    dene = zeros(length(lats), length(mod.alt));
    
    for li = 1:length(lats)
        dene(li, :) = interp_sami(mod, [lats(li), lon]);
    end
    nmf2(i, :) = max(dene, [], 2);
    for li = 1:length(lats)
        hmf2(i, li) = mean(mod.alt(dene(li, :) == nmf2(i, li)));
    end
    if i == 1
        dene_truth = dene;
        continue
    end
    pfreq_diff = (elec2freq(dene) - elec2freq(dene_truth)) / 1E3;

    
    %% plot
    if 0
    tiledlayout(1, 3)
    nexttile
    colormap parula
    hold on
    contourf(lats, mod.alt, elec2freq(dene_truth') / 1E3, 25)
    
    title('Truth')
    clim([0 6])
    colorbar
    xlim(latlimit)
    ylim(htlimit)
    xlabel('Lat (°)')
    ylabel('Alt (km)')
    colorbar
    for rxi = 1:2
        rays = homed_rays{1, rxi};
        for ri = 1:length(rays)
            for ri2 = 1:length(rays{ri})
                plot(rays{ri}(ri2).lat, rays{ri}(ri2).height, '-g', 'LineWidth', 2)
            end
        end
    end
    grid on
    grid minor
    hold off

    nexttile
    hold on
    contourf(lats, mod.alt, elec2freq(dene') / 1E3, 25)
    title('Observed')
    clim([0 6])
    colorbar
    xlim(latlimit)
    ylim(htlimit)
    xlabel('Lat (°)')
    % set(gca, 'YTickLabels', '')
    colorbar
    for rxi = 1:2
        rays = homed_rays{i, rxi};
        for ri = 1:length(rays)
            for ri2 = 1:length(rays{ri})
                plot(rays{ri}(ri2).lat, rays{ri}(ri2).height, '-m', 'LineWidth', 2)
            end
        end
    end
    grid on
    grid minor
    hold off

    nexttile
    hold on
    contourf(lats, mod.alt, pfreq_diff', 25)
    title('Difference')
    clim([-2 2])
    % set(gca, 'YTickLabels', '')
    xlim(latlimit)
        ylim(htlimit)

    xlabel('Lat (°)')
    hC = colorbar;
    ylabel(hC, 'Plasma Freq (MHz)', 'FontSize', 20)
    for rxi = 1:2
        rays = homed_rays{i, rxi};
        truth_rays = homed_rays{1, rxi};

        for ri = 1:length(truth_rays)
            for ri2 = 1:length(truth_rays{ri})
                plot(truth_rays{ri}(ri2).lat, truth_rays{ri}(ri2).height, '-g', 'LineWidth', 2)
            end
        end
        for ri = 1:length(rays)
            for ri2 = 1:length(rays{ri})
                plot(rays{ri}(ri2).lat, rays{ri}(ri2).height, '-m', 'LineWidth', 2)
            end
        end
    end
    grid on
    grid minor
    hold off
    end
end

hold off




%% HmF2 & NmF2
close
figure('units','normalized','outerposition',[0 0 1 1]);


tiledlayout(1, 2)
nexttile
hold on
for i = 1:length(modlist)
    plot(lats, elec2freq(nmf2(i, :))/1E3, 'LineWidth', 4)
end
hold off
ylim([0, 8])
grid on
grid minor
ylabel('foF2 (MHz)')
xlabel('Lat (°)')

nexttile
hold on
for i = 1:length(modlist)
plot(lats, hmf2(i, :), 'LineWidth', 4)
end
hold off
grid on
    grid minor
ylabel('hmF2 (km)')
xlabel('Lat (°)')
ylim([0, 400])
legend(modlist)













