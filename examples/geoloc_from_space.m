%% geoloc_from_space
% Identify the location of a surface transmitter from two or more
% measurements in space
% #1 Determine a search area on the surface
% #2 Home rays to the spacecraft from points in the search area
% #3 Compare against the observed time-difference of arrival


%% Inputs
in_fn_fmt = '/Users/chartat1/data/sami3/2017_tid/sami_mat/{YYYY-mm-dd_HHMM}.mat';
plot_dirn = '/Users/chartat1/data/itsi_sim/plots/';
geoloc_mov_plt = [plot_dirn, 'mov_geoloc/altlat_{YYYY-mm-dd-HHMM}.png'];
out_fn_fmt = '/Users/chartat1/data/itsi_sim/geoloc/sim_gs_rays_%s_%i.mat';

time = datenum(2017, 1, 10, 18, 0, 0);
aoa_times = time - 1/24:5/60/24:time + 1/24;

OX_mode = 1;
gdfreq = [10]; %, 15, 20, 100];
sounding_freqs = 2:0.1:8;

satlats = [35, 45];  % will move through these for Doppler plots
satlons = [282.5, 282.5] - 360;
satalt = 400;
gdlat_tx = 38;
gdlon = 282.5 -360;

maxdist = 1E5;  % meters from homing
tol = [1e-7 0.01 25];       % ODE solver tolerance and min max stepsizes
homing_tol_m = 750;

fprintf('Started a new raytracing expt\n')


%% Ground-to-satellite raytracing
sami = loadstruct(filename(in_fn_fmt, time));

% sami = subgrid_sami(sami, satlats, satlon);
txloc = [gdlat_tx, gdlon, 0];

[iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    Bx, By, Bz, geomag_grid_parms] = gen_grid_parms(sami);

for l = 1%:length(satlats)
    rxloc = [satlats(l), satlons(l), satalt];
    % raytrace
    homed_rays = gen_ionogram(gdfreq, OX_mode, txloc, rxloc, ...
        iono_en_grid, iono_en_grid_5, collision_freq, ...
        iono_grid_parms, geomag_grid_parms, Bx, By, Bz, ...
        maxdist, tol, homing_tol_m);
    savestruct(sprintf(out_fn_fmt, 'ground_space', l), homed_rays)
end


%% Calculate the differential delays
delays = zeros(length(satlats), length(gdfreq)) * NaN;

for l = 1:length(satlats)
    homed_rays = loadstruct(sprintf(out_fn_fmt, 'ground_space', l));
    for r = 1:length(homed_rays)
        delays(l, r) = homed_rays{r}.group_range_to_rx;
    end
end

obs_tdoa = (delays(2, :) - delays(1, :));


%% Perform a recursive grid search
recon_time = datenum(2017, 1, 10, 18, 5, 0);
recon_model = loadstruct(filename(in_fn_fmt, recon_time));

[iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    Bx, By, Bz, geomag_grid_parms] = gen_grid_parms(recon_model);
freq = 10;
searchalt = 0;
lat_delta = 1;
lon_delta = 0.5;
searchlats = [min(satlats) - 2.7:lat_delta:max(satlats) + 2.5];
searchlons = [min(satlons) - 1.7:lon_delta:max(satlons) + 1.5];
[slo2, sla2] = meshgrid(searchlons, searchlats);

npts = 15;
min_recon_tdoa_err = 1E6;
recon_tdoa_err_tol = 3;

while min_recon_tdoa_err > recon_tdoa_err_tol
    %% TDOA search
    recon_tdoa = tdoa_search(searchlats, searchlons, searchalt, ...
        satlats, satlons, satalt, npts, freq, OX_mode, ...
        iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
        Bx, By, Bz, geomag_grid_parms, maxdist, tol, homing_tol_m);
    
    %% Calculate errors
    recon_tdoa_err = abs(recon_tdoa - obs_tdoa);
    
    idx = recon_tdoa_err <= min(recon_tdoa_err(:) + recon_tdoa_err_tol);
    minlats = sla2(idx);
    minlons = slo2(idx);

    if min(recon_tdoa_err(:)) < recon_tdoa_err_tol
        min_recon_tdoa_err = min(recon_tdoa_err(:));
        [slo2, sla2] = meshgrid(searchlons, searchlats);
        fprintf('********** min(recon_tdoa_err): %1.1f km **********\n', min_recon_tdoa_err)
        fprintf('Search condition met\n')
        break
    end

    %% todo Account for multiple good-ish solutions
    searchlats = min(minlats) - lat_delta:lat_delta/10:max(minlats) + lat_delta;
    searchlons = min(minlons) - lon_delta:lon_delta/10:max(minlons) + lon_delta;
    [slo2, sla2] = meshgrid(searchlons, searchlats);
    lat_delta = (max(searchlats) - min(searchlats)) / 10;
    lon_delta = (max(searchlons) - min(searchlons)) / 10;

    assert(min_recon_tdoa_err > min(recon_tdoa_err(:)), 'No longer improving...')
    min_recon_tdoa_err = min(recon_tdoa_err(:));
    fprintf('********** min(recon_tdoa_err): %1.1f km **********\n', min_recon_tdoa_err)
end


%% Plot the cost function ellipse
close
subplot(2, 1, 1)
plot(searchlats, recon_tdoa_err(:, searchlons == gdlon), 'LineWidth', 3); 
xlabel("Lat (°)"); 
ylabel("TDOA err (km)");
grid on; grid minor
subplot(2, 1, 2)
plot(searchlons, recon_tdoa_err(searchlats == gdlat_tx, :), 'LineWidth', 3); 
xlabel("Lon (°)"); 
ylabel("TDOA err (km)");
grid on; grid minor

figure;
hold on
contourf(searchlons, searchlats, recon_tdoa_err)
plot(gdlon, gdlat_tx, '.g', 'MarkerSize', 30)
plot(minlon, minlat, '.m', 'MarkerSize', 30)
plot(satlons, satlats, '.r', 'MarkerSize', 30)
plot(slo2(:), sla2(:), '.k')
xlabel('Lon (°)')
ylabel('Lat (°)')
legend({'TDOA Err. (km)', 'True Location', 'Closest guess', 'Spacecraft locations', 'Search points'})
xlim([min(searchlons), max(searchlons)])
ylim([min(searchlats), max(searchlats)])
cb = colorbar;
ylabel(cb,'TDOA err (km)','FontSize',16,'Rotation',90)
hold off

dist_km = distance([minlat, minlon], [gdlat_tx, gdlon], wgs84Ellipsoid) / 1E3;
target_tdoa_err = recon_tdoa_err(sla2 == gdlat_tx & slo2 == gdlon);
fprintf('Min TDOA err: %1.1f km\n', min(recon_tdoa_err(:)))
fprintf('TDOA err at the target: %1.1f km\n', target_tdoa_err)

fprintf('Distance to target: %1.1f km\n', dist_km)

% TODO: print TDOA error at transmitter location

%% Plot the ionosphere at the two times
lati = sami.lat >= min(satlats) & sami.lat <= max(satlats);
loni = sami.lon == gdlon + 360;
lats = sami.lat(lati);
lons = sami.lon(loni);
ne = zeros(length(sami.alt), sum(lati), 3);
ne(:, :, 1) = squeeze(sami.dene(:, lati, loni));
ne(:, :, 2) = squeeze(recon_model.dene(:, lati, loni));
ne(:, :, 3) = ne(:, :, 2) - ne(:, :, 1);
for ct = 1:3
    subplot(1, 3, ct)

    contourf(lats, sami.alt, squeeze(ne(:, :, ct)))
    if ct < 3
    clim([0 4E5])
    else
        clim([-5E4 5E4])
    end
    colorbar
    xlabel('Lat (°)')
    ylabel('Alt (km)')
end



%% Define TDOA search algorithm
function recon_tdoa = tdoa_search(searchlats, searchlons, searchalt, ...
    satlats, satlons, satalt, npts, freq, OX_mode, ...
    iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    Bx, By, Bz, geomag_grid_parms, maxdist, tol, homing_tol_m)

% Calculate differential TDOAs at grid search points
verbose = false;
recon_tdoa = zeros(length(searchlats), length(searchlons)) * NaN;
for l1 = 1:length(searchlats)
    
    for l2 = 1:length(searchlons)
        fprintf('l1: %i of %i, l2: %i of %i\n', l1, length(searchlats), l2, length(searchlons))
        txloc = [searchlats(l1), searchlons(l2), searchalt];
        delays = zeros(length(satlats)) * NaN;
        for l = 1:length(satlats)
            rxloc = [satlats(l), satlons(l), satalt];
            % raytrace
            homed_rays = gen_ionogram(freq, OX_mode, txloc, rxloc, ...
                iono_en_grid, iono_en_grid_5, collision_freq, ...
                iono_grid_parms, geomag_grid_parms, Bx, By, Bz, ...
                maxdist, tol, homing_tol_m, verbose);
            min_delay = 1E9;
            for r = 1:length(homed_rays)
                if homed_rays{1}(r).group_range_to_rx < min_delay
                    min_delay = homed_rays{1}(r).group_range_to_rx;
                end
            end
            delays(l) = min_delay;

        end
        delays(delays == 1E9) = NaN;
        recon_tdoa(l1, l2) = delays(2) - delays(1);
    end
end


end

























