%% ground_space_sim.m
% Simulate transionospheric sounding for Nebula sponsor

%% Inputs
time = datenum(2017, 1, 10, 18, 0, 0);
in_fn_fmt = '/Users/chartat1/data/sami3/2017_tid/sami_mat/{YYYY-mm-dd_HHMM}.mat';

OX_mode = 0;
gdfreq = [10];
aoafreq = 5;
freqs = 2:0.1:9;
satlats = 30:0.1:55;  % will move through these for Doppler plots
satlon = 282.5 - 360;
satalt = 580;
gdlat = 38;
gdlon = 282.5 - 360;

maxdist = 1E5;  % meters from homing
tol = [1e-7 0.01 25];       % ODE solver tolerance and min max stepsizes
homing_tol_m = 1500;

fprintf('Started a new raytracing expt\n')


%% Ground-to-satellite raytracing
sami = loadstruct(filename(in_fn_fmt, time));
sami_5 = loadstruct(filename(in_fn_fmt, time + 5/60/24));
sami = subgrid_sami(sami, satlats, satlon);
sami_5 = subgrid_sami(sami_5, satlats, satlon);
txloc = [gdlat, gdlon, 0];


% First do the case where the receiver is exactly overhead
rxloc = [gdlat, gdlon, satalt];
homed_rays = gen_ionogram(freqs, OX_mode, txloc, rxloc, ...
    sami.iono_en_grid, sami_5.iono_en_grid, sami.collision_freq, ...
    sami.iono_grid_parms, sami.geomag_grid_parms, ...
    sami.Bx, sami.By, sami.Bz, ...
    maxdist, tol, homing_tol_m, true, false);

lon = satlon;

%%
close all
figure('units','normalized','outerposition',[0 0 1 1]);
colormap parula
hold on
[~, hC] = contourf(sami.lat, sami.alt, elec2freq(sami.dene(:, :, sami.lon == lon)) / 1E3, 50);
set(hC, 'LineStyle', 'none')
xlabel('Lat (°)')
ylabel('Alt (km)')
xlim([37, 39])
ylim([0, 600])
cl = colorbar;
ylabel(cl, 'Plasma Frequency (MHz)')
set(gca, 'color', 'k', 'FontSize', 30, 'XColor', 'w', 'YColor', 'w')
set(cl, 'XColor', 'w', 'YColor', 'w')
set(gcf, 'InvertHardCopy', 'off');
set(gcf,'Color',[0 0 0]); % RGB values [0 0 0] indicates black color

for ri = 1:length(homed_rays)
    for r = 1:length(homed_rays{ri})
        plot(homed_rays{ri}(r).lat, homed_rays{ri}(r).height, '-w', 'LineWidth', 2)
    end
end
plot(txloc(1), txloc(3), '.g', 'MarkerSize', 50)
plot(rxloc(1), rxloc(3), '.r', 'MarkerSize', 50)

%% 
close all
figure

fof2 = squeeze(elec2freq(max(sami.dene, [], 1)) / 1E3);
contourf(sami.lon, sami.lat, fof2, 50)
hold on

for ri = 1%:length(homed_rays)
    for r = 1:length(homed_rays{ri})
        plot(homed_rays{ri}(r).lon, homed_rays{ri}(r).lat, '-w', 'LineWidth', 2)
    end
end
xlim([-78, -77])
ylim([37, 39])
fprintf('FoF2 at transmitter: %1.1f \n', elec2freq(max(sami.dene(:, sami.lat == 38, sami.lon == gdlon))) / 1E3)


%% Then do a more realistic overpass, with 1 s per frequency
ct = 0;
for l = 1:length(freqs)
    rxloc = [satlats(l), gdlon, satalt];
    % raytrace
    homed_rays{l} = gen_ionogram(freqs(l), OX_mode, txloc, rxloc, ...
        sami.iono_en_grid, sami_5.iono_en_grid, sami.collision_freq, ...
        sami.iono_grid_parms, sami.geomag_grid_parms, ...
        sami.Bx, sami.By, sami.Bz, ...
        maxdist, tol, homing_tol_m);
end

%% 
close all
figure('units','normalized','outerposition',[0 0 1 1]);
colormap parula
hold on
[~, hC] = contourf(sami.lat, sami.alt, elec2freq(sami.dene(:, :, sami.lon == lon)) / 1E3, 50);
set(hC, 'LineStyle', 'none')
xlabel('Lat (°)')
ylabel('Alt (km)')
xlim([34, 39])
ylim([0, 600])
cl = colorbar;
ylabel(cl, 'Plasma Frequency (MHz)')
set(gca, 'color', 'k', 'FontSize', 30, 'XColor', 'w', 'YColor', 'w')
set(cl, 'XColor', 'w', 'YColor', 'w')
set(gcf, 'InvertHardCopy', 'off');
set(gcf,'Color',[0 0 0]); % RGB values [0 0 0] indicates black color

for ri = 1:length(homed_rays)
    for r = 1:length(homed_rays{ri})
        for r2 = 1:length(homed_rays{ri}{r})
        plot(homed_rays{ri}{r}(r2).lat, homed_rays{ri}{r}(r2).height, '-w', 'LineWidth', 2)
        plot(homed_rays{ri}{r}(r2).rxloc(1), homed_rays{ri}{r}(r2).rxloc(3), '.r', 'MarkerSize', 50)
        end

    end
end
plot(txloc(1), txloc(3), '.g', 'MarkerSize', 50)



%% 
figure 
hold on
colormap cool

for ri = 1:length(homed_rays)
    for r = 1:length(homed_rays{ri})
        for r2 = 1:length(homed_rays{ri}{r})

            freq = homed_rays{ri}{r}(r2).frequency;
            rg = homed_rays{ri}{r}(r2).group_range_to_rx;
            elevation = 0;


            scatter(freq, rg, 50, elevation, 'filled', 'o')
        end

    end
end
xlim([2, 9])
ylim([0, 1500])
set(gca, 'YDir','reverse')
ylabel({'Virtual Range (km)'})
xlabel('Tx Freq (MHz)')
set(gcf, 'InvertHardCopy', 'off');

grid on
grid minor
