function D = load_wspr_csv(csv_fname, stime, etime, fmin, fmax, snrmin, distmin, distmax)
%% Load the WSPR data from https://www.wsprnet.org/drupal/downloads
% table:
% 1 spot ID, 2 Timestamp, 3 Reporter, 4 Reporter grid, 5 SNR, 
% 6 Frequency (MHz), 7 Tx callsign, 8 Tx grid, 9 Tx power (dBm), 
% 10 Drift (Hz/min), 11 Great circle distance (km), 12 Azimuth (deg), 
% 13 Band (-1: LF, 0: MF, 1: 160m, 3: 80m, 5: 60m, 7: 40m, 10: 30m, ...).
% 14 Version
% 15 Code (0 = good)
% 
% Test inputs
% csv_fname = '~/Downloads/wsprspots-2014-05.csv';
% stime = datetime(2014, 5, 23, 6, 0, 0);
% etime = datetime(2014, 5, 23, 12, 0, 0);
% fmin = 2;
% fmax = 30;
% snrmin = -25;
% distmin = 500;
% distmax = 2000;

%% Load
M = readtable(csv_fname);

%% read/convert
clear D; 
Flags = table2array(M(:, 15));
SNR = table2array(M(:, 5));

D.times = table2array(M(:, 2)) / 86400 + datenum(1970, 1, 1);
D.txloc = table2array(M(:, 8));
D.rxloc = table2array(M(:, 4));
D.freq = table2array(M(:, 6));
D.dist = table2array(M(:, 11));
txpower = table2array(M(:, 9));


% Select relevant vals
idx = (Flags == 0) & (D.times >= datenum(stime)) & (D.times <= datenum(etime)) & ...
    (D.freq >= fmin) & (D.freq <= fmax) & (SNR >= snrmin) & ...
    (D.dist >= distmin) & (D.dist <= distmax) & (cellfun(@length, D.txloc) == 6) ...
    & (cellfun(@length, D.rxloc) == 6);

assert(~isempty(idx), 'No valid data - change parameters')

fn = fieldnames(D);
for k = 1:numel(fn)
    D.(fn{k}) = D.(fn{k})(idx);
end

% Convert from Maidenhead to lat/lon
for t = 1:length(D.times)
    [D.txlat(t), D.txlon(t)] = maidenhead_to_ll(D.txloc{t});
    [D.rxlat(t), D.rxlon(t)] = maidenhead_to_ll(D.rxloc{t});
end

D = rmfield(D, {'txloc', 'rxloc'});


























