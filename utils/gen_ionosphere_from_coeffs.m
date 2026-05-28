function [iono_en_grid, iono_grid_parms] = gen_ionosphere_from_coeffs( ...
    X, alts, lats, lons)
%% gen_ionosphere_from_coeffs
% X = 3 * len(lats) * len(lons), 1st set are nmf2, then hmf2, then H
% X = 1 is 1E6 el. cm-3 nmf2, 300 km hmf2, 100 km H (scale height)

%% def iono_grid_parms
lat_step = unique(diff(lats));
lon_step = unique(diff(lons));
alt_step = unique(diff(alts));

iono_grid_parms = [...
    min(lats), lat_step, length(lats), ...
    min(lons), lon_step, length(lons), ...
    min(alts), alt_step, length(alts)];

%% Generate ionosphere from coefficients
shape = [length(lats), length(lons)];
npts = length(lats) * length(lons);

Np = reshape(X(1:npts), shape);  % peak density coefficients
hp = reshape(X(npts + 1:2*npts), shape); % peak height coefficients
H = reshape(X(2*npts + 1:end), shape);  % topside scale height coefficients
c = 1;
iono_en_grid = calc_chapman_profs(lats, lons, alts, Np * 1E6, hp * 300, H * 100, c);
