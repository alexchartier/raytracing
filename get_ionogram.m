function [ionogram, idx] = get_ionogram(homed_rays, sl_cutoff_pct)
%% get_ionogram

%% Check if we're vertical or oblique
if sum(homed_rays{1}(1).txloc == homed_rays{1}(1).rxloc) == 3
    vert = true;
else
    vert = false;
end

%% Grab the freqs and ranges
freqs = [];
ranges = [];
for ri = 1:length(homed_rays)
    for r = 1:length(homed_rays{ri})
        freqs = [freqs, homed_rays{ri}(r).frequency];
        ranges = [ranges, homed_rays{ri}(r).group_range_to_rx];
    end
end

% Pull out the unique ones (rounded to nearest km)
[ranges, ids] = unique(round(ranges));
freqs = freqs(ids);

% discard anything inside near-straightline cutoff range
sep = dist_txrx(homed_rays{1}(1).txloc, homed_rays{1}(1).rxloc);
rg_cutoff = sep * (1 + sl_cutoff_pct / 100);
rgi = ranges > rg_cutoff;
ranges = ranges(rgi);
freqs = freqs(rgi);
[ranges, si] = sort(ranges);
freqs = freqs(si);


%% Set output params
ionogram.freqs = freqs;
ionogram.ranges = ranges;
ionogram.mode = homed_rays{1}(1).OX_mode;
ionogram.N_local_tx = homed_rays{1}(1).electron_density(1);
ionogram.N_local_rx = interp1(homed_rays{1}(1).geometric_distance, ...
    homed_rays{1}(1).electron_density, ...
    homed_rays{1}(1).geometric_dist_to_rx);

ionogram.N_local = (ionogram.N_local_tx + ionogram.N_local_rx) / 2;
ionogram.sep = dist_txrx(homed_rays{1}(1).txloc, homed_rays{1}(1).rxloc);

if ionogram.sep == 0
    ionogram.local_gradient = 0;
else
    ionogram.local_gradient = (ionogram.N_local_rx - ionogram.N_local_tx) / ionogram.sep;
end

ionogram.rg_cutoff = rg_cutoff;
ionogram.txloc = homed_rays{1}(1).txloc; 
ionogram.rxloc = homed_rays{1}(1).rxloc;
[ionogram.lat_ref, ionogram.lon_ref] = midpointLatLon(...
    ionogram.txloc(1), ionogram.txloc(2), ionogram.rxloc(1), ionogram.rxloc(2));
 
end






