%% vis_synth.m
% Generate vertical ionograms

%% Set inputs
clear

ray_fn = '~/Downloads/vert_rays.mat';
txloc = [50.2, 0.5, 0];
time = datenum(2017, 1, 10, 18, 0, 0);
OX_mode = [-1, 0, 1];
freqs = 2:0.1:15;

% pharlap parameters
maxdist = 1E5;  % meters from homing
tol = [1e-7 0.01 25];       % ODE solver tolerance and min max stepsizes
homing_tol_m = 500;

R12 = 150;
UT = [2017, 1, 10, 18, 0];
alts = 90:3:1000;

%% Generate ionosphere
[iono, iono_extra] = iri2020(txloc(1), txloc(2), R12, UT, alts(1), ...
    alts(2) - alts(1), length(alts));
Ne = iono(1, :) / 1E6;

model.time = time;
model.alt = alts;
model.lat = [txloc(1) - 5, txloc(1) + 5];
model.lon = [txloc(2) - 5, txloc(2) + 5];
model.dene = repmat(Ne', [1, 2, 2]);


[iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    Bx, By, Bz, geomag_grid_parms] = gen_grid_parms(model);

%% Run raytracing
for oxi = 1:length(OX_mode)
    homed_rays{oxi} = gen_ionogram(freqs, OX_mode(oxi), txloc, txloc, ...
        iono_en_grid, iono_en_grid_5, collision_freq, ...
        iono_grid_parms, geomag_grid_parms, Bx, By, Bz, ...
        maxdist, tol, homing_tol_m);
end

homed_rays{1}{1}.iono_model = model;
savestruct(ray_fn, homed_rays)

%% Plot
homed_rays = loadstruct(ray_fn);
close


subplot(1, 2, 1)
colors = {'r', 'm', 'b'};
hold on
for oxi = 1:length(OX_mode)
    for ri = 1:length(homed_rays{oxi})
        for r = 1:length(homed_rays{oxi}{ri})
            plot(homed_rays{oxi}{ri}(r).frequency, ...
                homed_rays{oxi}{ri}(r).group_range_to_rx / 2, ...
                sprintf('.%s', colors{oxi}), 'MarkerSize', 10)
        end
    end
end
text(6,900, 'X', 'Color', 'r', 'FontSize', 14); 
text(6,850, 'No B', 'Color', 'm', 'FontSize', 14); 
text(6,800, 'O', 'Color', 'b', 'FontSize', 14); 
text(6,950, 'Iono plasmafreq', 'Color', 'k', 'FontSize', 14); 
plot(elec2freq(Ne) / 1E3, alts, 'k', 'LineWidth', 2)
xlabel('Freq (MHz)')
ylabel('Virtual Height (km)')
ylim([0, 1000])
grid on
grid minor

subplot(1, 2, 2)
plot(homed_rays{2}{end}.lat, homed_rays{2}{end}.height)
xlabel('Lat (°)')
ylabel('Height (km)')
xlim([txloc(1) - 2, txloc(1) + 2])
ylim([0, 1000])
grid on
grid minor
