function homed_rays = gen_ionogram(freqs, OX_mode, txloc, rxloc, ...
    iono_en_grid, iono_en_grid_5, collision_freq, ...
    iono_grid_parms, geomag_grid_parms, Bx, By, Bz, ...
    maxdist, tol, homing_tol_m, varargin)

%% generate an ionogram using raytrace_itsi
% % Note: assumes Bx/By/Bz are same size as iono_en_grid
%
% % Usage:
% homed_rays = gen_ionogram(freqs, OX_mode, txloc, rxloc, ...
%     iono_en_grid, iono_en_grid_5, collision_freq, ...
%     iono_grid_parms, geomag_grid_parms, Bx, By, Bz, ...
%     maxdist, tol, homing_tol_m, varargin);
%
% TODO: add the time to the output

%% Verbose switch
if length(varargin) >= 1
    verbose = varargin{1};
else
    verbose = true;
end

if length(varargin) == 2
    straightline_break = varargin{2};
else
    straightline_break = true;
end

%% -180/180 the lons
if txloc(2) > 180
    txloc(2) = txloc(2) - 360;
end
if rxloc(2) > 180
   rxloc(2) = rxloc(2) - 360;
end


%% subgrid the ionosphere for speed
lats_full = arange(iono_grid_parms(1), iono_grid_parms(2), iono_grid_parms(3));
lons_full = arange(iono_grid_parms(4), iono_grid_parms(5), iono_grid_parms(6));

lati = lats_full >= min([txloc(1), rxloc(1)]) - iono_grid_parms(2) * 2 & ...
    lats_full <= max([txloc(1), rxloc(1)]) + iono_grid_parms(2) * 2;
loni = lons_full >= min([txloc(2), rxloc(2)]) - iono_grid_parms(5) * 2 & ...
    lons_full <= max([txloc(2), rxloc(2)]) + iono_grid_parms(5) * 2;

assert(sum(lati) > 0, 'Transmitter or receiver outside grid')
assert(sum(loni) > 0, 'Transmitter or receiver outside grid')
lats = lats_full(lati);
lons = lons_full(loni);

iono_en_grid = iono_en_grid(lati, loni, :);
iono_en_grid_5 = iono_en_grid_5(lati, loni, :);
collision_freq = collision_freq(lati, loni, :);
Bx = Bx(lati, loni, :);
By = By(lati, loni, :);
Bz = Bz(lati, loni, :);

iono_grid_parms = [ ...
    min(lats), iono_grid_parms(2), length(lats), ...
    min(lons), iono_grid_parms(5), length(lons), ...
    iono_grid_parms(7), iono_grid_parms(8), iono_grid_parms(9)];
geomag_grid_parms = [ ...
    min(lats), geomag_grid_parms(2), length(lats), ...
    min(lons), geomag_grid_parms(5), length(lons), ...
    geomag_grid_parms(7), geomag_grid_parms(8), geomag_grid_parms(9)];

fof2_max = elec2freq(max(max(max(iono_en_grid / 1E6)))); 

%% sanity checks
assert(txloc(1) >= min(lats) && txloc(1) <= max(lats) && ...
    txloc(2) >= min(lons) && txloc(2) <= max(lons), ...
    'Transmitter outside ionosphere grid')
assert(rxloc(1) >= min(lats) && rxloc(1) <= max(lats) && ...
    rxloc(2) >= min(lons) && rxloc(2) <= max(lons), ...
    'Receiver outside ionosphere grid')

%% guesstimate the likely raypath geometry
% az/el/range from transmitter to receiver
[gc_az, gc_el, sep] = geodetic2aer( ...
    rxloc(1), rxloc(2), rxloc(3), ...
    txloc(1), txloc(2), txloc(3), ...
    wgs84Ellipsoid("km"));

% same for transmitter to midpoint
[latMid, lonMid] = midpointLatLon(txloc(1), txloc(2), rxloc(1), rxloc(2));
[~, tangent_el, ~] = geodetic2aer( ...
    latMid, lonMid, 0, ...
    txloc(1), txloc(2), txloc(3), ...
    wgs84Ellipsoid("km"));

[~, f_el, ~] = geodetic2aer( ...
    latMid, lonMid, 300, ...
    txloc(1), txloc(2), txloc(3), ...
    wgs84Ellipsoid("km"));


if verbose
    fprintf('\n')
    fprintf('%1.1f °N, %1.1f °E, %1.0f km\n', txloc(1), txloc(2), txloc(3))
    fprintf('Tangent elevation: %1.1f, 300-km tangent: %1.1f\n', tangent_el, f_el)
    fprintf('Separation: %1.1f km\n', sep)
end
% ad-hoc classification of NVIS vs oblique
%TODO: do a better job of selecting az/el for the ground-to-space case
if sep < 100 % NVIS
    azarr = 0:10:360; %-10:5:10; %0:5:360
else % oblique
    fac = ceil(30*sqrt(100/sep));
    azarr = linspace(gc_az-fac, gc_az+fac, 20);
    azarr(azarr > 360) = azarr(azarr > 360) - 360;
end

if txloc(3) > 300 && rxloc(3) > 300 % topside sounding case
    elvarr = linspace(tangent_el, gc_el, 100);

elseif (txloc(3) < 100 && rxloc(3) > 100) || ...
        (txloc(3) > 100 && rxloc(3) < 100) % transionospheric case
    elvarr = gc_el - 20:gc_el + 20;
elseif (txloc(3) < 100 && rxloc(3) < 100) % ground case
    elvarr = 0:90;
end

elvarr = elvarr(elvarr >= -90);

if txloc(3) == 0
    elvarr = elvarr(elvarr > 0);
end

%% Loop over freqs
% cut the loop if we've got to a quasi-straight line ray
ct = 1;
homed_rays = {};
if verbose
    fprintf('foF2: %1.1f MHz\n', fof2_max)
end
for f = 1:length(freqs)
    if verbose
        fprintf('f: %1.1f MHz\n', freqs(f))
    end

    %fprintf("Min el.: %1.1f Max el.: %1.1f\n", min(elvarr(:)), max(elvarr(:)))
    ray = ...
        raytrace_itsi(freqs(f), OX_mode, txloc, rxloc, ...
        iono_en_grid, iono_en_grid_5, collision_freq, ...
        iono_grid_parms, Bx, By, Bz, geomag_grid_parms, ...
        elvarr, azarr, maxdist, tol, 0, homing_tol_m, verbose);

    if isstruct(ray)
        homed_rays{ct} = ray;
        ct = ct + 1;

        longest_ray_dist = 0;
        for ri = 1:length(ray)
            if longest_ray_dist < ray(ri).geometric_dist_to_rx
                longest_ray_dist = ray(ri).geometric_dist_to_rx;
            end
        end

        % break out if we're really close to straight
        if (longest_ray_dist / sep) < 1.005 && straightline_break 
            fprintf('Breaking due to straight line prop.: %1.1f MHz\n', freqs(f))
            return
        end
    end

    % break out if we're vertical sounding above fof2
    if (sep == 0) && (freqs(f) > fof2_max + 1) && ~isstruct(ray) 
        fprintf('Breaking due to >1MHz above max foF2 in vertical sounding mode: %1.1f MHz, %1.1f MHz, \n', freqs(f), fof2_max)
        return
    end
end


%% check that at least one ray went through the ionosphere
goodone = 0;
for r = 1:length(homed_rays)
    for ri = 1:length(homed_rays{r})
        if sum(homed_rays{r}(ri).electron_density) > 0
            goodone = 1;
        end
    end
end
if ~goodone && verbose
    fprintf('No good rays found\n')
end






















