

%% Set inputs
sami_fn = 'sami3_regulargrid_elec_density_2014May23.nc';
out_fn = 'sim_ois_rays.mat';
time = datenum(2014, 5, 23, 1, 0, 0);
alts = 90:2:450;
OX_mode = 0;
maxdist = 5E3;
tol = [1e-7 0.01 25];
refractive_ind = 1;
nhops = 1;
freqs = 2:0.1:20;

txloc = [62.3, -145.3, 0];  % NOME
% txloc = [-23.66, 144.14, 0];  % JORN 
rxloc = [64.5, -165.4, 0];


%% Load ionosphere from SAMI3
sami_vars = {'alt', 'lat', 'lon', 'dene0', 'time'};
sami = load_sami(sami_fn, sami_vars);

%% loop over times
disp(time)
rays = [];

%% Get ionosphere at time t from SAMI3 into PHaRLAP format
sami_t = interp_sami(sami, alts, datetime(time, 'ConvertFrom', 'datenum'));
[iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    Bx, By, Bz, geomag_grid_parms] = gen_grid_parms(sami_t);

%% Load links at time t
nlinks = length(freqs);

%% Raytrace
rays = [];
for l = 1:nlinks
    [ray, elvarr, azarr] = raytrace_3dhome(freqs(l), txloc, rxloc, ...
        iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
        Bx, By, Bz, geomag_grid_parms, OX_mode, tol, refractive_ind, nhops, maxdist);
    rays(l).ray = ray;
    rays(l).home = ray.home;
    rays(l).txloc = txloc;
    rays(l).rxloc = rxloc;
end

%% Plotting rays on globe
figure
hold on
satglobe4e
% iono_panel = surf(squeeze(X), squeeze(Y), squeeze(Z), 'FaceColor', 'none', 'EdgeColor', 'none');
% alpha = 0.7; % globe transparency level, 1 = opaque, through 0 = invisible
% set(iono_panel, 'FaceColor', 'texturemap', 'CData', squeeze(V), 'FaceAlpha', alpha, 'EdgeColor', 'none');

% plot the rays
Re = 6371E3;
for r = 1:length(rays)
    ray = rays(r).ray;
    if ~ray.home
        continue
    end

    if ~isempty(ray.initial_elev)

        sph = [ray.height * 1E3 + Re; deg2rad(ray.lat); deg2rad(ray.lon)];
        cart = sphcart(sph');
        plot3(cart(:, 1), cart(:, 2), cart(:, 3), 'w')
    end
end


hold off

%%
%
% set(gca,'color','black')
% export_fig('~/Downloads/twosat_sim.tif');


%% Plotting 2 - ionogram plot
figure
colormap('cool')

xlimit = [2, 20];
ylimit = [500, 1500];
colorlimit = [-130 -80];
ct = 0;
hold on
for r = 1:length(rays)
    ray = rays(r).ray;
    if ray.home
        freq = ray.frequency;
        rg = ray.group_range(end);
        loss = ray.absorption(end) + fspl(rg * 1E3, 3E8 / (freq * 1E6));
        scatter(freq, rg, 150, -loss, 'filled', 'o')
    end
end

ylim(ylimit)
xlim(xlimit)
ylabel({'Virtual Range (km)'})
    xlabel('Tx Freq (MHz)')

%legend({'X', 'No B', 'O'})
% title(sprintf('foF2: %1.1f MHz, hmF2: %1.1f km, S/C alt: %1.1f km', D.fof2(timeidx), D.hmf2(timeidx), alt_iss))


hold off

grid on
grid minor




























