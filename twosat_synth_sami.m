%% twosat_synth_sami.m
% Simulate two satellite topside sounder mission data using the SAMI3 model
%
%NOTE: segfault sometimes happens if you launch too soon after startup. 
% Could be somehow related to parallel pool startup?
% Run twosat_synth.m first for a minute or two to avoid this. 


clear

%% inputs
in_fn_fmt = 'data/sami/{YYYY-mm-dd_HHMM}.mat';
out_fn_fmt = 'data/sami_osse/sim_twosat_rays_%s_%i.mat';

time = datenum(2015, 3, 4, 8, 0, 0);
freqs = 2:0.1:10; % 2:0.2:20; % 20;
OX_modes = -1; 

satlats = 40:80;  % will move through these
sep = 5;
satlon = 172.5;
satlon(satlon > 180) = satlon(satlon > 180) - 360;
satalt = 600;

% homing
elvarr = -70:2:0;
azarr = -10:5:10;
maxdist = 1E5;  % meters from homing
blind_range = 30E3; % two-way (e.g. there-and-back range)
tol = [1e-7 0.01 25];       % ODE solver tolerance and min max stepsizes
homing_tol_m = 750; 

fprintf('Started a new raytracing expt\n')


%% Load ionosphere
% NOTE: not sure there's any point in reducing the grid sizing in terms of
% speed
sami = loadstruct(filename(in_fn_fmt, time));

lati_1 = sami.lat >= min(satlats);
lati_2 = sami.lat <= max(satlats);
loni = find(sami.lon == satlon);
alti = 55:2:length(sami.alt);  % cut off at 200km
sami.alt = sami.alt(alti);
sami.lat = sami.lat(lati_1 & lati_2);
sami.lon = sami.lon(loni - 1:loni + 1);
sami.dene = sami.dene(alti, lati_1 & lati_2, loni - 1:loni + 1);
[iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    Bx, By, Bz, geomag_grid_parms] = gen_grid_parms(sami);
lats = arange(iono_grid_parms(1), iono_grid_parms(2), iono_grid_parms(3));
lons = arange(iono_grid_parms(4), iono_grid_parms(5), iono_grid_parms(6));
alts = arange(iono_grid_parms(7), iono_grid_parms(8), iono_grid_parms(9));

%% loop over latitude
for lati = 1:range(satlats) - sep
    %% Specify transmitter and receiver locs
    txloc = [satlats(lati), satlon, satalt];
    rxloc = [satlats(lati + 5), satlon, satalt];

    %% Calculate foF2 for reference
    iono_en_prof = nan(size(alts));
    for i = 1:length(alts)
        iono_en_prof(i) = interp2(lats, lons, ...
            squeeze(iono_en_grid(:, :, alts == alts(i)))',...
            txloc(1), txloc(2));
    end
    fprintf('FoF2: %1.2f\n', sqrt(80.6 * max(iono_en_prof)./ 1E6))

    %% Loop over O/X mode
    for oxi = 1:length(OX_modes)
        OX_mode = OX_modes(oxi);
        switch OX_mode
            case -1
                OX_mode_name = 'X';
            case 1
                OX_mode_name = 'O';

            case 0
                OX_mode_name = 'No_B';
        end

        out_fn = sprintf(out_fn_fmt, OX_mode_name, lati);

        %% Produce an ionogram through raytracing
        homed_rays = gen_ionogram(freqs, OX_mode, txloc, rxloc, ...
            iono_en_grid, iono_en_grid_5, collision_freq, ...
            iono_grid_parms, Bx, By, Bz, geomag_grid_parms, ...
            elvarr, azarr, maxdist, tol, homing_tol_m);

        %% save
        clear raytrace_3d
        homed_rays{1}(1).iono_en_grid = iono_en_grid;
        homed_rays{1}(1).iono_lat = lats;
        homed_rays{1}(1).iono_lon = lons;
        homed_rays{1}(1).iono_alt = alts;

        savestruct(out_fn, homed_rays)
        fprintf("saved %i homed_rays to %s\n", length(homed_rays), out_fn)

    end
end

% for i = 1:length(homed_ray)
%     disp(homed_ray(i).group_range_to_rx - tr(i).group_range_to_rx)
% end

% %% Convert to netCDF (not working yet)
% % loop over lat
% for lati = 1:range(satlats) - sep
%     % Loop over O/X mode
%     for OX_mode = -1:2:1
%         in_fn = sprintf(out_fn_fmt, OX_mode_name, lati);
%         nc_fn = sprintf(nc_fn_fmt, OX_mode_name, lati);
% 
%         homed_rays = loadstruct(in_fn);
%         struct2nc(homed_rays, nc_fn);
%     end
% end

% 
%% Plot rays on electron density grid
clf
hold on
% freq
% contourf(homed_rays{1}(1).iono_lat, homed_rays{1}(1).iono_alt, ...
%     squeeze(elec2freq(homed_rays{1}(1).iono_en_grid(:, 2, :)))'/1E3, 50); 
contourf(homed_rays{1}(1).iono_lat, homed_rays{1}(1).iono_alt, ...
    squeeze(homed_rays{1}(1).iono_en_grid(:, 2, :))', 50); 

fmin = sqrt(min(freqs));
fmax = sqrt(max(freqs));
fd = fmax - fmin; 

for lati = 1:range(satlats) - sep

    homed_rays = loadstruct(sprintf(out_fn_fmt, OX_mode_name, lati));

    for i = 1:length(homed_rays)
        for j = 1:length(homed_rays{i})
            ray = homed_rays{i}(j);
            fac = (sqrt(ray.frequency)  - fmin) / fd;
            hi = ray.height <= satalt;
            plot(ray.lat(hi), ray.height(hi), 'color', [1, fac, fac], 'LineWidth', 1.2)
%             plot(ray.lat(1), ray.height(1), '.m', 'markersize', 20)
%             plot(ray.rxloc(1), ray.rxloc(3), '.m', 'markersize', 20)
        end
    end
end
xlim([45, 73]); 
clim(gca, [0, 3E5]);
h = colorbar;
ylabel(h, 'Electron Density (el. cm^{-3})')
xlabel('Latitude (degrees)')
ylabel('Altitude (km)')
set(gca, 'FontSize', 36)


%% Plotting 2 - ionogram plot
for lati =3 % 1:range(satlats) - sep
    clf
    homed_rays = loadstruct(sprintf(out_fn_fmt, OX_mode_name, lati));
    plot_ionogram(homed_rays)
    pause(1)
end
