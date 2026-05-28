%% hiamcm_raytracing
% Raytrace through TID ionosphere

clear
%% set inputs
in_fn_fmt = '/Users/chartat1/data/sami3/2017_tid/sami_mat/{YYYY-mm-dd_HHMM}.mat';
plot_dirn = '/Users/chartat1/data/sami3/2017_tid/plots/';
geoloc_mov_plt = [plot_dirn, 'mov_geoloc_fliptxrx/altlat_{YYYY-mm-dd-HHMM}.png'];
out_fn_fmt = '/Users/chartat1/data/sami3/2017_tid/recon_rays_ground_space/O_mode/sim_twosat_rays_%s_%i.mat';
out_fn_fmt_vert_iono = '/Users/chartat1/data/sami3/2017_tid/recon_rays/O_mode/vert_{YYYY-mm-dd_HHMM}_%1.1fN_%1.1fE_%ikm.mat';
out_table_fn = '/Users/chartat1/data/sami3/2017_tid/ranges.xls';
out_fn_fmt_oblique_iono = ['/Users/chartat1/data/sami3/2017_tid/recon_rays/', ...
    'O_mode/oblique_{YYYY-mm-dd_HHMM}_%1.1fN_%1.1fE_%1.1fN_%1.1fE_%ikm.mat'];

time = datenum(2017, 1, 10, 18, 0, 0);
aoa_times = time - 1/24:5/60/24:time + 1/24;

OX_mode = 1;
gdfreq = [10];
aoafreq = 5;
sounding_freqs = 2:0.1:15;

satlats = 30:55;  % will move through these for Doppler plots
sounding_lats = 30:55; % 34.1:46.1;  % will move through these for topside sounding
satlon = 282.5 - 360;
satalt = 580;
alongtrack_spacing_km = 600;
gdlat_tx = 38;
gdlon = 282.5 - 360;
gdlat_rx = 39;

maxdist = 1E5;  % meters from homing
tol = [1e-7 0.01 25];       % ODE solver tolerance and min max stepsizes
homing_tol_m = 1500;

fprintf('Started a new raytracing expt\n')


%% Ground-to-ground raytracing for AoA demo
az_out = zeros(length(aoa_times), 5) * NaN;
el_out = zeros(length(aoa_times), 5) * NaN;
txloc = [gdlat_tx, gdlon, 0];
rxloc = [gdlat_rx, gdlon, 0];

for t = 1:length(aoa_times)
    %%
    sami = loadstruct(filename(in_fn_fmt, aoa_times(t)));

    aoa_rays = gen_ionogram(aoafreq, OX_mode, txloc, rxloc, ...
        sami.iono_en_grid, sami.iono_en_grid, sami.collision_freq, ...
        sami.iono_grid_parms, sami.geomag_grid_parms, ...
        sami.Bx, sami.By, sami.Bz, ...
        maxdist, tol, homing_tol_m);
    for r = 1:length(aoa_rays{1})
        az_out(t, r) = aoa_rays{1}(r).initial_bearing;
        el_out(t, r) = aoa_rays{1}(r).initial_elev;
    end

    %% plotting
    close
    figure('Position', [1200, 800, 1200, 800])
    latidx = sami.lat >= 30 & sami.lat <=50;
    dene = squeeze(sami.dene(:, latidx, ismember(sami.lon - 360, gdlon)));
    hold on
    plot(rxloc(1), rxloc(3), '.g', 'MarkerSize', 50)
    plot(txloc(1), txloc(3), '.r', 'MarkerSize', 50)
    [~, hC] = contourf(sami.lat(latidx), sami.alt, elec2freq(dene) / 1E3);
    for r = 1:length(aoa_rays{1})
        plot(aoa_rays{1}(r).lat, aoa_rays{1}(r).height, '-m', 'LineWidth', 3)
    end
    hold off
    legend({'Transmitter', 'Receiver'})
    xlim([37, 40])
    ylim([0, 400])
    xlabel('Lat (°)')
    ylabel('Alt (km)')
    % set(hC, 'LineStyle', 'none')
    clim([0, 6])
    a = colorbar;
    a.Label.String = 'Plasma Freq (MHz)';
    title(sprintf('%s, %1.1f MHz, %1.1f° E', ...
        filename('{YYYY-mm-dd HH:MM UT}', aoa_times(t)), aoafreq, gdlon))

    export_fig(filename(geoloc_mov_plt, aoa_times(t)))
    pause(0.1)

end

% more plotting
az_out(az_out == 0) = NaN;
el_out(el_out == 0) = NaN;
az_out_2 = az_out;
az_out_2(el_out > 90) = az_out_2(el_out > 90) + 180;
clf
subplot(2, 1, 1)
hold on
for i = 1:5
    plot(aoa_times, az_out_2(:, i), '-k', 'LineWidth', 3)
end
ylabel('Received Azimuth (°)')
grid on
grid minor
hold off

set(gca, 'XTickLabels', '')
subplot(2, 1, 2)
hold on
for i = 1:5
    plot(aoa_times, el_out(:, 1), '-k', 'LineWidth', 3)
end
hold off
ylabel('Received Elevation (°)')
datetick('keeplimits')
grid on
grid minor
xlabel('Time (UT)')


%% satellite topside sounding (vertical)
sami = loadstruct(filename(in_fn_fmt, time));

for l = 1:length(sounding_lats)
    txloc = [sounding_lats(l), satlon, satalt];
    rxloc = txloc; 
    homed_rays = gen_ionogram(sounding_freqs, OX_mode, txloc, rxloc, ...
        sami.iono_en_grid, sami.iono_en_grid, sami.collision_freq, ...
        sami.iono_grid_parms, sami.geomag_grid_parms, ...
        sami.Bx, sami.By, sami.Bz, ...
        maxdist, tol, homing_tol_m);
    out_fn = sprintf(filename(out_fn_fmt_vert_iono, time), txloc(1), txloc(2), txloc(3));
    savestruct(out_fn, homed_rays)
    fprintf('Saved to %s\n', out_fn)
end


%% satellite topside sounding (oblique)
sami = loadstruct(filename(in_fn_fmt, time));

for l = 1:length(sounding_lats) - 1
    txloc = [sounding_lats(l), satlon, satalt];
    refloc = [sounding_lats(l + 1), satlon, satalt];
    rxloc = calc_rx_loc_for_ois(txloc, refloc, alongtrack_spacing_km);

    homed_rays = gen_ionogram(sounding_freqs, OX_mode, txloc, rxloc, ...
        sami.iono_en_grid, sami.iono_en_grid, sami.collision_freq, ...
        sami.iono_grid_parms, sami.geomag_grid_parms, ...
        sami.Bx, sami.By, sami.Bz, ...
        maxdist, tol, homing_tol_m);

    out_fn = sprintf(filename(out_fn_fmt_oblique_iono, time), ...
        txloc(1), txloc(2), rxloc(1), rxloc(2), txloc(3));
    savestruct(out_fn, homed_rays)
    fprintf('Saved to %s\n', out_fn)
end


%% Ionospheric reconstruction
sami = loadstruct(filename(in_fn_fmt, time));
homed_rays = loadstruct(sprintf(out_fn_fmt, 'vert', 2));
recon_vert_ionogram(homed_rays, sami)


%% Ground-to-satellite raytracing
sami = loadstruct(filename(in_fn_fmt, time));
sami_5 = loadstruct(filename(in_fn_fmt, time + 5/60/24));
sami = subgrid_sami(sami, satlats, satlon);
sami_5 = subgrid_sami(sami_5, satlats, satlon);
txloc = [gdlat_tx, gdlon, 0];


for l = 1:length(satlats)
    rxloc = [satlats(l), gdlon, satalt];
    % raytrace    
    homed_rays = gen_ionogram(gdfreq, OX_mode, txloc, rxloc, ...
        sami.iono_en_grid, sami_5.iono_en_grid, sami.collision_freq, ...
        sami.iono_grid_parms, sami.geomag_grid_parms, ...
        sami.Bx, sami.By, sami.Bz, ...
        maxdist, tol, homing_tol_m);
    savestruct(sprintf(out_fn_fmt, 'ground_space', l), homed_rays)
end


%% Plotting the ground-to-space raytracing
ranges = zeros(length(gdfreq), length(satlats), 4) * NaN;
dopplers = zeros(length(gdfreq), length(satlats), 4) * NaN;
seps = zeros(length(gdfreq), 1);
sami = loadstruct(filename(in_fn_fmt, time));

figure('Position', [1200, 800, 1200, 800])
latidx = sami.lat >= min(satlats) & sami.lat <= max(satlats);
dene = squeeze(sami.dene(:, latidx, ismember(sami.lon - 360, gdlon)));
hold on
[~, hC] = contourf(sami.lat(latidx), sami.alt, elec2freq(dene) / 1E3, 50);
set(hC, 'LineStyle', 'None')
[~, hC] = contourf(sami.lat(latidx), [0, min(sami.alt)], zeros(2, sum(latidx)), 50);
set(hC, 'LineStyle', 'None')

% store out range and Doppler, and generate the ray plot
for l = 1:length(satlats)
    homed_rays = loadstruct(sprintf(out_fn_fmt, 'ground_space', l));
    txloc = homed_rays{1}.txloc; 
    rxloc = homed_rays{1}.rxloc; 
    [~, ~, seps(l)] = geodetic2aer( ...
        rxloc(1), rxloc(2), rxloc(3), ...
        txloc(1), txloc(2), txloc(3), ...
        wgs84Ellipsoid("km"));
    for r = 1:length(homed_rays)
        ray = homed_rays{r};
        for r1 = 1:length(ray)
            if r1 > size(ranges, 3)
                continue
            end
            ranges(gdfreq == ray(r1).frequency, l, r1) = ray(r1).group_range_to_rx;
            dopplers(gdfreq == ray(r1).frequency, l, r1) = ray(r1).final_Doppler;
        end
    end

    % ray plot
        for r = 1:length(homed_rays)
            ray = homed_rays{r};
        for r1 = 1:length(ray)

        plot(ray(r1).lat, ray(r1).height, '-w', 'LineWidth', 3)
        end
        end

    plot(rxloc(1), rxloc(3), '.r', 'MarkerSize', 10)
    plot(txloc(1), txloc(3), '.g', 'MarkerSize', 10)
end
xlim([min(satlats), max(satlats)])
ylim([0, 280])
    a = colorbar;
    a.Label.String = 'Plasma Freq (MHz)';
    xlabel('Lat (°)')
    ylabel('Alt (km)')
hold off


%% plot the range-Doppler
nnan = sum(sum(isnan(ranges), 3), 2 );
maxnan = size(ranges, 2) * size(ranges, 3);
goodfreqs = gdfreq(nnan < maxnan);
goodrg = ranges(nnan < maxnan, :, :);
gooddop = dopplers(nnan < maxnan, :, :);
figure
colormap jet
for f = 2:length(goodfreqs)
    subplot(length(goodfreqs), 1, f)
    hold on
    for mode = 1:size(goodrg, 3)
        scatter(satlats, goodrg(f, :, mode) - seps', 40, gooddop(f, :, mode), 'filled')
    end
    ylabel(sprintf('%i MHz\nIono Rg. (km)', goodfreqs(f)))
    if f == length(goodfreqs)
        xlabel('S/C lat (°)')
    else
        set(gca, 'XTicklabels', '')
    end
    grid on
    grid minor
        a = colorbar;
    a.Label.String = 'Doppler (Hz)';
    hold off
end

rg_10MHz = ranges(2, :, 1);
rg_15MHz = ranges(3, :, 1);
rg_20MHz = ranges(4, :, 1);
rg_100MHz = ranges(5, :, 1);
% T = table({'Sat Lat' '10MHz', '15MHz', '20MHz', '100 MHz'}, satlats', rg_10MHz', rg_15MHz', rg_20MHz', rg_100MHz');
% writetable(T, out_table_fn)









