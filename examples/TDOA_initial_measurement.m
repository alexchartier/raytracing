%% inputs 
time = datenum(2017, 1, 11, 11, 0, 0);

sc_crd_fn1 = '~/data/nebula/STK_positions/tdoa_sats//CubeSat1_Fixed_Position_Velocity.txt';
sc_crd_fn2 = '~/data/nebula/STK_positions/tdoa_sats//CubeSat2_Fixed_Position_Velocity.txt';
gd_loc = [38, -77.5, 0];

run_name = 'vert';

switch run_name
    case {'truth'}
        ionosphere_map_fn_fmt = ...
            '/Users/chartat1/data/sami3/2017_tid/sami_mat/{YYYY-mm-dd_HHMM}.mat';
    case {'vert'}
        ionosphere_map_fn_fmt = ...
            '/Users/chartat1/data/sami3/2017_tid/recon_iono/vert_recon_{YYYY-mm-dd_HHMM}*.mat';
    case {'oblique'}
        ionosphere_map_fn_fmt = ...
            '/Users/chartat1/data/sami3/2017_tid/recon_iono/oblique_recon_{YYYY-mm-dd_HHMM}*.mat';
end
tdoa_file_fmt = '~/data/nebula/TDOA_test/%s_TDOA_Profile_%iMHz.xlsx';
ray_fn_fmt = '~/data/nebula/TDOA_test/rays_%s_%iMHz.mat';
Re = 6371E3;
freq = 10;

%% Read
ionosphere_map_fname = glob(filename(ionosphere_map_fn_fmt, time));
ionosphere_map_fname = ionosphere_map_fname{1};
CubeSat1_file = readmatrix(sc_crd_fn1);
CubeSat2_file = readmatrix(sc_crd_fn2);
ground_loc_x_y_z = sphcart([gd_loc(3) + Re, deg2rad(gd_loc(1)), deg2rad(gd_loc(2))]); % x,y,z (m)

%% Calculate positions
CubeSat1.x_pos = CubeSat1_file(:,5) * 1000; % m
CubeSat1.y_pos = CubeSat1_file(:,6) * 1000; % m
CubeSat1.z_pos = CubeSat1_file(:,7) * 1000; % m
receiver1_loc_x_y_z = [CubeSat1.x_pos, CubeSat1.y_pos, CubeSat1.z_pos]; % m
% receiver1_loc_all = cartsph(receiver1_loc_x_y_z);
% receiver1_loc_all(:, 1) = receiver1_loc_all(:, 1) - Re;
% receiver1_loc_all(:, 2) = rad2deg(receiver1_loc_all(:, 2));
% receiver1_loc_all(:, 3) = rad2deg(receiver1_loc_all(:, 3));

% CubeSat2
CubeSat2.x_pos = CubeSat2_file(:,5) * 1000; % m
CubeSat2.y_pos = CubeSat2_file(:,6) * 1000; % m
CubeSat2.z_pos = CubeSat2_file(:,7) * 1000; % m
receiver2_loc_x_y_z = [CubeSat2.x_pos, CubeSat2.y_pos, CubeSat2.z_pos]; % m

% Calculate delta T from ionosphere (Alex code) (us)

TDOA = [];
n = [];
rays = {};
for i = 540:40:660
    [ionosphere_delay_us1, rays_1] = calc_ionospheric_delay(...
        ionosphere_map_fname, ground_loc_x_y_z, receiver1_loc_x_y_z(i,:), freq);
    [ionosphere_delay_us2, rays_2] = calc_ionospheric_delay(...
        ionosphere_map_fname, ground_loc_x_y_z, receiver2_loc_x_y_z(i,:), freq);
    
    rays.sc1{i} = rays_1;
    rays.sc2{i} = rays_2;


    if (isnan(ionosphere_delay_us1(1)))
        ionosphere_delay_us1 = 0;
    end
    if (isnan(ionosphere_delay_us2(1)))
        ionosphere_delay_us2 = 0;
    end
    if (isnan(ionosphere_delay_us1(1)) && isnan(ionosphere_delay_us2(1) == 0))
        n = [n 1e6]; % us
        TDOA = [TDOA (ionosphere_delay_us1(1) - ionosphere_delay_us2(1))];
    else
        n = [n 1e-3]; % us
        TDOA = [TDOA (ionosphere_delay_us1(1) - ionosphere_delay_us2(1))];
    end
end

% rxloc_1 = cartsph(receiver1_loc_x_y_z(i, :));
% rxloc_1 = [rad2deg(rxloc_1(2)), rad2deg(rxloc_1(3)), rxloc_1(1)/1E3 - 6371];
% rxloc_2 = cartsph(receiver2_loc_x_y_z(i, :));
% rxloc_2 = [rad2deg(rxloc_2(2)), rad2deg(rxloc_2(3)), rxloc_2(1)/1E3 - 6371];
% earth_example; plot_rays(rays_1, rxloc_1, rxloc_2); plot_rays(rays_2, rxloc_1, rxloc_2)


data = [540:40:660; TDOA; n];

writematrix(data',sprintf(tdoa_file_fmt, run_name, freq))
savestruct(sprintf(ray_fn_fmt, run_name, freq), rays)



%% Plotting/analysis
truth_iono = '/Users/chartat1/data/sami3/2017_tid/sami_mat/{YYYY-mm-dd_HHMM}.mat';
recon_iono = glob(filename(...
    '/Users/chartat1/data/sami3/2017_tid/recon_iono/vert_recon_{YYYY-mm-dd_HHMM}*.mat', ...
    time));
ray_fn_fmt = '~/data/nebula/TDOA_test/rays_%s_%iMHz.mat';
truth_rays = loadstruct(sprintf(ray_fn_fmt, 'truth', freq));
recon_rays = loadstruct(sprintf(ray_fn_fmt, 'vert', freq));

D1 = loadstruct(filename(truth_iono, time));
D2 = loadstruct(recon_iono{1});
D1lon = D1.lon;
D1lon(D1lon > 180) = D1lon(D1lon > 180) - 360;
[~, ~, id1] = closest(D1lon, gd_loc(2));
[~, ~, id2] = closest(D2.lon, gd_loc(2));
truth_slice = elec2freq(squeeze(D1.dene(:, :, id1))) / 1E3;
recon_slice = elec2freq(squeeze(D2.dene(:, :, id2))) / 1E3;
 
truth_nmf2 = squeeze(max(D1.dene, [], 1));
recon_nmf2 = squeeze(max(D2.dene, [], 1));
close
climit = [min(recon_nmf2(:)), max(recon_nmf2(:))];
subplot(2, 1, 1)
hold on
contourf(D1.lon, D1.lat, truth_nmf2, 350)
ylim([30, 45])
xlim([-80, -75] + 360)
plot(gd_loc(2) + 360, gd_loc(1), 'xr', 'MarkerSize', 50)
clim(climit)
colorbar

subplot(2, 1, 2)
hold on
contourf(D2.lon, D2.lat, recon_nmf2, 30)
ylim([30, 45])
xlim([-80, -75])
clim(climit)
colorbar

plot(gd_loc(2), gd_loc(1), 'xr', 'MarkerSize', 50)

%%
close

climit = [0; max(recon_slice(:))];
subplot(2, 1, 1)
[~, hc] = contourf(D1.lat, D1.alt, truth_slice, 50);
set(hc, 'LineStyle', 'none')
xlim([30, 45])
clim(climit)

hold on
for r = 540:40:660
    rays = truth_rays.sc1{r};
    for ri = 1:length(rays)
        plot(rays{ri}.lat, rays{ri}.height, '-w')
    end

    rays = recon_rays.sc1{r};
    for ri = 1:length(rays)
        plot(rays{ri}.lat, rays{ri}.height, '-m')
    end
end
colorbar

subplot(2, 1, 2)
hold on
[~, hc] = contourf(D2.lat, D2.alt, recon_slice, 50);
set(hc, 'LineStyle', 'none')
for r = 540:40:660
    rays = recon_rays.sc1{r};
    for ri = 1:length(rays)
        plot(rays{ri}.lat, rays{ri}.height, '-w')
    end

end
xlim([30, 45])
clim(climit)
colorbar























