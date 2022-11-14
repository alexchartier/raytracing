function raytrace_sami(year, month, day, hours, sami_fn, wspr_fn, out_fn_fmt, varargin)
% Function to raytrace through the SAMI3 model and compare against WSPR 
% reported links. Monthy WSPR csv files are available here: 
% http://www.wsprnet.org/drupal/downloads
% 
% % Example inputs
% % Required
% year = 2014;
% month = 5;
% day = 23; 
% hours = 0:23;
% sami_fn = 'sami3_regulargrid_elec_density_2014May23.nc';
% wspr_fn = 'wsprspots-2014-05.csv';
% out_fn_fmt = 'links_{yyyymmmdd-HHMM}.nc';
% 
% % Optional (w. defaults):
% % Grid altitudes
% alts = 90:2:450;
% 
% % Prop details
% OX_mode = 0;
% maxdist = 5E3;  % Tolerance for homing (meters)
% tol = [1e-7 0.01 25];       % ODE solver tolerance and min max stepsizes (see PHaRLAP docs)
% refractive_ind = 1;  % at the transmitter location - set to 1 if transmitting outside iononsphere
% nhops = 1;
% 
% % Analysis parameters
% fmin = 2;
% fmax = 30;
% snrmin = -25;  
% % Note this is not a true "SNR" - seems to be ratio of
% % narrow carrier band to a broader "test band". Therefore it is usually negative but
% % higher values may indicate "radio quiet" at the station. 
% distmin = 500;
% distmax = 2000;
% npts_wanted = 200; % links at each time
% 
% % Call
% raytrace_sami(year, month, day, hours, sami_fn, wspr_fn, out_fn_fmt)
% raytrace_sami(year, month, day, hours, sami_fn, wspr_fn, out_fn_fmt, ...
%   'fmax', 20, 'snrmin', -30)


%% Set default args
params = {'alts', 'OX_mode', 'maxdist', 'tol', 'refractive_ind', 'nhops' ...
    'fmin', 'fmax', 'snrmin', 'distmin', 'distmax', 'npts_wanted'};
defaults = {90:2:450, 0, 5E3, [1e-7 0.01 25], 1, 1, 2, 30, -25, 500, 2000, 200};
varargparse(varargin, params, defaults);


%% Set the analysis times
times = datetime(year, month, day, hours, 0, 0);


%% Load ionosphere from SAMI3 and monthly links from WSPR
sami = load_sami(sami_fn);
wspr = load_wspr_csv(wspr_fn, min(times), max(times), fmin, fmax, ...
    snrmin, distmin, distmax);


%% loop over times
for t = 1:length(times)
    time = times(t);
    disp(time)
    links = [];

    %% Get ionosphere at time t from SAMI3 into PHaRLAP format
    sami_t = interp_sami(sami, alts, time);
    [iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
        Bx, By, Bz, geomag_grid_parms] = gen_grid_parms(sami_t);


    %% Get WSPR data close to the right time
    wspr_ti = abs(wspr.times - datenum(time)) < 5/60/24;
    fn = fieldnames(wspr);
    wspr_t = [];
    for k = 1:numel(fn)
        wspr_t.(fn{k}) = wspr.(fn{k})(wspr_ti);
    end

    wspr_out = declump_wspr(wspr_t, npts_wanted);

    txlocs = [wspr_out.txlat; wspr_out.txlon; ones([1, length(wspr_out.times)])]';
    rxlocs = [wspr_out.rxlat; wspr_out.rxlon; ones([1, length(wspr_out.times)])]';
    freqs = wspr_out.freq;


    %% Load links at time t
    nlinks = length(txlocs);
    links.txlocs = txlocs;
    links.rxlocs = rxlocs;
    links.freqs = freqs;
    links.home = zeros([nlinks, 1]);
    links.ray = cell([nlinks, 1]);
        
        
    %% Raytrace
    for l = 1:nlinks
        [ray, elvarr, azarr] = raytrace_3dhome(freqs(l), txlocs(l, :), rxlocs(l, :), ...
            iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
            Bx, By, Bz, geomag_grid_parms, OX_mode, tol, refractive_ind, nhops, maxdist);
        links.ray{l} = ray;
        links.home(l) = ray.home;
    end

    struct2nc(links, filename(out_fn_fmt, times(t)), 'netcdf4_classic');
end




% %% Plotting
% rotang = -30;
% 
% txloc = links{2}.txlocs;
% rxloc = links{2}.rxlocs;
% home = links{2}.home;
% 
% latI_bad = [];
% lonI_bad = [];
% latI_good = [];
% lonI_good = [];
% for l = 1:length(home)
%     % Great circle calcs
%     [latI,lonI] = interpm([txloc(l, 1), rxloc(l, 1)], ...
%         [txloc(l, 2), rxloc(l, 2)], 1, 'gc');
%     if home(l) == 0
%         latI_bad = [latI_bad; latI; NaN];
%         lonI_bad = [lonI_bad; lonI; NaN];
%     else
%         latI_good = [latI_good; latI; NaN];
%         lonI_good = [lonI_good; lonI; NaN];
%     end
% end
% 
% rotang = -30;
% earth_example
% % Plot the links
% CartL = sphcart([64E5 * ones(size(latI_bad)), deg2rad(latI_bad), deg2rad(lonI_bad)]);
% h3 = plot3(CartL(:, 1), CartL(:, 2), CartL(:, 3), '-r', 'LineWidth', 3);
% rotate(h3, [0 0 1], rotang);
% 
% 
% CartL = sphcart([64E5 * ones(size(latI_good)), deg2rad(latI_good), deg2rad(lonI_good)]);
% h3 = plot3(CartL(:, 1), CartL(:, 2), CartL(:, 3), '-g', 'LineWidth', 3);
% rotate(h3, [0 0 1], rotang);
% %CartL = sphcart([64E5 * ones(size(latI_f)), deg2rad(latI_f), deg2rad(lonI_f)]);
% 
% rotate(globe, [0 0 1], rotang)


























