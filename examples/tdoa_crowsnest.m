%% inputs 
clear
time = datenum(2017, 1, 11, 11, 5, 0);
% 
% sat_locs = [
%      36, -80.5, 250; 
%     37, -75.5, 250; 
%     40, -80.5, 250; 
%     41, -75.5, 250; 
% ];


sat_locs = [
    31, -77.5, 200; 
    38, -77.5, 200; 
];

sat_locs_high = [
    28, -77.5, 400; 
    43, -77.5, 400; 
];

txloc = [35, -77.5, 0];
ionosphere_map_fn_fmt = ...
            '/Users/chartat1/data/sami3/2017_tid/sami_mat/{YYYY-mm-dd_HHMM}.mat';
ray_fn_fmt = '~/data/crowsnest/TDOA_test/rays_%s_%iMHz.mat';
Re = 6371E3;
freq = 4;

maxdist = 1E5;  % meters from homing
tol = [1e-7 0.01 25];       % ODE solver tolerance and min max stepsizes
homing_tol_m = 750;
OX_modes = [-1, 1];
c = 299792458;

%% Read
ionosphere_map_fname = glob(filename(ionosphere_map_fn_fmt, time));
ionosphere_map_fname = ionosphere_map_fname{1};
model = loadstruct(ionosphere_map_fname);


%% Raytrace
for i = 1:size(sat_locs, 1)
    rxloc = sat_locs(i, :);
    for mi = 1:length(OX_modes)

        rays = gen_ionogram(freq, OX_modes(mi), txloc, rxloc, ...
            model.iono_en_grid, model.iono_en_grid, model.collision_freq, ...
            model.iono_grid_parms, model.geomag_grid_parms, ...
            model.Bx, model.By, model.Bz, ...
            maxdist, tol, homing_tol_m);
        if length(rays) > 0
            low_rays{i, mi} = rays{1};
        end
    end
end


for i = 1:size(sat_locs_high, 1)
    rxloc = sat_locs_high(i, :);
    for mi = 1:length(OX_modes)

        rays = gen_ionogram(freq, OX_modes(mi), txloc, rxloc, ...
            model.iono_en_grid, model.iono_en_grid, model.collision_freq, ...
            model.iono_grid_parms, model.geomag_grid_parms, ...
            model.Bx, model.By, model.Bz, ...
            maxdist, tol, homing_tol_m);
        if length(rays) > 0
            high_rays{i, mi} = rays{1};
        end
    end
end


%% Raytrace for 450 km case



%% Plotting/analysis

modlon = model.lon;
modlon(modlon > 180) = modlon(modlon > 180) - 360;
[~, ~, id1] = closest(modlon, txloc(2));
truth_slice = elec2freq(squeeze(model.dene(:, :, id1))) / 1E3;


close

climit = [0; 3];
% subplot(2, 1, 1)
[~, hc] = contourf(model.lat, model.alt, truth_slice, 50);
set(hc, 'LineStyle', 'none')
xlim([25, 45])
ylim([0, 500])
xlabel('Lat (°)')
ylabel('Alt (km)')
clim(climit)

hold on
for r = 1:size(low_rays, 1)
    for ri = 1:size(low_rays, 2)
        try
            size(low_rays{r, ri});
            if ri == 1
                plot(low_rays{r, ri}.lat, low_rays{r, ri}.height, '-m')
            else
                plot(low_rays{r, ri}.lat, low_rays{r, ri}.height, '--m')
            end
            plot(sat_locs(r, 1), sat_locs(r, 3), '.g', 'MarkerSize', 50)
        catch 'MS'
        end
    end

end

for r = 1:size(high_rays, 1)
    for ri = 1:size(high_rays, 2)
        try
            size(high_rays{r, ri});
            if ri == 1
                plot(high_rays{r, ri}.lat, high_rays{r, ri}.height, '-m')
            else
                plot(high_rays{r, ri}.lat, high_rays{r, ri}.height, '--m')
            end
            plot(sat_locs_high(r, 1), sat_locs_high(r, 3), '.r', 'MarkerSize', 50)
        catch 'MS'
        end
    end

end

h = colorbar;
ylabel(h, 'Plasma Freq (MHz)')
% 
% subplot(2, 1, 2)
% hold on
% [~, hc] = contourf(D2.lat, D2.alt, recon_slice, 50);
% set(hc, 'LineStyle', 'none')
% for r = 540:40:660
%     rays = recon_rays.sc1{r};
%     for ri = 1:length(rays)
%         plot(rays{ri}.lat, rays{ri}.height, '-w')
%     end
% 
% end
% xlim([30, 45])
% clim(climit)
% colorbar























