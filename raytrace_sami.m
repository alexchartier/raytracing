%% raytrace_sami.m
% Script to raytrace through the SAMI3 model and compare against WSPR 
% reported links. Monthy WSPR csv files are available here: 
% http://www.wsprnet.org/drupal/downloads
% 

%% inputs
times = datetime(2014, 5, 23, 23:24, 0, 0);
sami_fn = 'sami3_regulargrid_elec_density_2014May23.nc';
wspr_fn = 'wsprspots-2014-05.csv';

% Grid altitudes
alts = 90:2:450;

% Prop details
OX_mode = 0;
maxdist = 5E3;  % meters from homing
tol = [1e-7 0.01 25];       % ODE solver tolerance and min max stepsizes
refractive_ind = 1;  % at the transmitter location
nhops = 1;

% Analysis parameters
fmin = 2;
fmax = 30;
snrmin = -25;
distmin = 500;
distmax = 2000;


%% Load ionosphere from SAMI3 and monthly links from WSPR
sami = load_sami(sami_fn);
wspr = load_wspr_csv(csv_fname, min(times), max(times), fmin, fmax, ...
    snrmin, distmin, distmax);


%% loop over times
links = cell([length(times), 1]);
for t = 1:length(times)
    time = times(t);

    %% Get ionosphere at time t from SAMI3 into PHaRLAP format
    sami_t = interp_sami(sami, alts, time);
    [iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
        Bx, By, Bz, geomag_grid_parms] = gen_iono_geomag_grids(sami_t);
    iono_pf_grid = sqrt(80.6 * iono_en_grid ./ 1E6);


    %% Get WSPR data close to the right time
    wspr_ti = abs(wspr.times - datenum(time)) < 5/60/24;
    fn = fieldnames(wspr);
    wspr_t = [];
    for k = 1:numel(fn)
        wspr_t.(fn{k}) = wspr.(fn{k})(wspr_ti);
    end

    txlocs = [wspr_t.txlat; wspr_t.txlon; ones([1, length(wspr_t.times)])]';
    rxlocs = [wspr_t.rxlat; wspr_t.rxlon; ones([1, length(wspr_t.times)])]';
    freqs = wspr_t.freq;


    %% Load links at time t
    nlinks = length(txlocs);
    links{t}.txlocs = txlocs;
    links{t}.rxlocs = rxlocs;
    links{t}.freqs = freqs;
    links{t}.home = zeros([nlinks, 1]);
    links{t}.ray = cell([nlinks, 1]);
        
        
    %% Raytrace
    for l = 1:nlinks
        [ray, elvarr, azarr] = raytrace_3dhome(freqs(l), txlocs(l, :), rxlocs(l, :), ...
            iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
            Bx, By, Bz, geomag_grid_parms, OX_mode, tol, refractive_ind, nhops, maxdist);
        links{t}.ray{l} = ray;
        links{t}.home(l) = ray.home;
    end
end


%% Plotting
rotang = -30;

txloc = links{2}.txlocs;
rxloc = links{2}.rxlocs;
home = links{2}.home;

latI_bad = [];
lonI_bad = [];
latI_good = [];
lonI_good = [];
for l = 1:length(home)
    % Great circle calcs
    [latI,lonI] = interpm([txloc(l, 1), rxloc(l, 1)], ...
        [txloc(l, 2), rxloc(l, 2)], 1, 'gc');
    if home(l) == 0
        latI_bad = [latI_bad; latI; NaN];
        lonI_bad = [lonI_bad; lonI; NaN];
    else
        latI_good = [latI_good; latI; NaN];
        lonI_good = [lonI_good; lonI; NaN];
    end
end

rotang = -30;
earth_example
% Plot the links
CartL = sphcart([64E5 * ones(size(latI_bad)), deg2rad(latI_bad), deg2rad(lonI_bad)]);
h3 = plot3(CartL(:, 1), CartL(:, 2), CartL(:, 3), '-r', 'LineWidth', 3);
rotate(h3, [0 0 1], rotang);


CartL = sphcart([64E5 * ones(size(latI_good)), deg2rad(latI_good), deg2rad(lonI_good)]);
h3 = plot3(CartL(:, 1), CartL(:, 2), CartL(:, 3), '-g', 'LineWidth', 3);
rotate(h3, [0 0 1], rotang);
%CartL = sphcart([64E5 * ones(size(latI_f)), deg2rad(latI_f), deg2rad(lonI_f)]);

rotate(globe, [0 0 1], rotang)


























