%% raytrace_sami.m
% Script to raytrace through the SAMI3 model and compare against WSPR 
% reported links. Monthy WSPR csv files are available here: 
% http://www.wsprnet.org/drupal/downloads
% 

%% inputs
times = datenum(2019, 2, 28):1/24:datenum(2019, 6, 1);

% Grid 
alts = 90:2:450;
lats = -90:-72;
lons = 151:180;
lons(lons > 180) = lons(lons > 180) - 360;

% Prop details
OX_mode = 0;
maxdist = 1E3;  % meters homing error (max)
tol = [1e-7 0.01 25];  % ODE solver tolerance and min max stepsizes
refractive_ind = 1;  % at the transmitter location (leave as 1 unless transmitting from within a plasma)
nhops_max = 2;
freqs = [4.1, 5.1 6.0 6.4 7.2]; 

% locations
txloc = [-77.8, 166.4, 0];
rxloc = [-89.9, 166.4, 1];

% Sunspot number: TODO get a real R12 set from the time
R12 = 100;


%% loop over times
links = cell([length(times), 1]);
for t = 1:length(times)
    time = times(t);
    UT = (time - floor(time)) * 24;

    %% Load ionosphere
    [iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
        Bx, By, Bz, geomag_grid_parms] = ...
        gen_iono_geomag_grids(alts, lats, lons, time, R12);

    %% Set up storage
    nlinks = length(freqs);
    links{t}.txloc = txloc;
    links{t}.rxloc = rxloc;
    links{t}.freqs = freqs;
    links{t}.home = zeros([nlinks, 1]);
        
        
    %% Raytrace
    for l = 1:nlinks
        [ray, elvarr, azarr] = raytrace_3dhome(freqs(l), txloc, rxloc, ...
            iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
            Bx, By, Bz, geomag_grid_parms, OX_mode, tol, refractive_ind,...
            nhops_max, maxdist);
        links{t}.ray(l) = ray;
        links{t}.home(l) = ray.home;
    end
end


%% Plotting
plot_rays(txloc, rxloc, links{1}.ray, [0 1 0], 180);























