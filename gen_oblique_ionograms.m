%% gen_oblique_ionograms.m
% Generate an input set of oblique ionograms with 600-km spacing 

%% inputs
sami_fn_fmt = '/Users/chartat1/data/sami3/2017_tid/sami_mat/{YYYY-mm-dd_HHMM}.mat';
sc_fn_fmt = ['~/data/nebula/STK_positions/Nebula_Constellation_LLA_AllTimesteps/', ...
    'Nebula_Constellation_LLA_Time%i.csv'];

out_fn_fmt = ['/Users/chartat1/data/sami3/gs_ionograms_3d/hiamcm/', ...
    '{YYYY-mm-dd_HHMM}/oblique/ionogram_%1.1fN_%1.1fE_to_%1.1fN_%1.1fE_%ikm.mat'];

times = datenum(2017, 1, 10, 12, 0, 0):3/24:datenum(2017, 1, 13);
OX_mode = 1;
freqs = 2:0.1:25;

alongtrack_spacing_km = 600;  % spacing between spacecraft

% pharlap parameters
maxdist = 1E5;  % meters from homing
tol = [1e-7 0.01 25];       % ODE solver tolerance and min max stepsizes
homing_tol_m = 1500;
fprintf('Started a new raytracing expt\n')


%% Calculate the spacecraft positions
nt = 67;
spacing = 6;
poslist = zeros(nt, 3);
az = zeros(nt - 1, 1);

wgs84 = wgs84Ellipsoid;
for t = 1:nt
    data = readtable(sprintf(sc_fn_fmt, t));
    XYZ  = lla2ecef([data.Var2, data.Var3, data.Var4 * 1E3], 'WGS84');

    if t == 1
        idx = find(data.Var2 > 20 & data.Var2 < 25);
        % idx = find(data.Var2 > -80 & data.Var2 < -68);
        idx = idx(2);
    end
    poslist(t, :) = XYZ(idx, :); % just grab the first one above 20°
end

SphV = cartsph(poslist);
sc_rad = mean(SphV(:, 1));
alt = sc_rad /1E3 - 6371;
arclen = 360 * alongtrack_spacing_km * 1E3 / (2 * pi * sc_rad);
txlats = rad2deg(SphV(:, 2));
txlons = rad2deg(SphV(:, 3));
rxlats = zeros(size(az));
rxlons = zeros(size(az));

for t = 1:nt-1
    az(t) = azimuth(txlats(t), txlons(t), txlats(t + 1), txlons(t + 1), wgs84, 'degrees');
    [rxlats(t), rxlons(t)] = reckon(txlats(t), txlons(t), arclen, az(t));
end
txlats = txlats(1:end-1);
txlons = txlons(1:end-1);

txlats = txlats(1:spacing:end);
txlons = txlons(1:spacing:end);
rxlats = rxlats(1:spacing:end);
rxlons = rxlons(1:spacing:end);


%% repmat out to cover all longitudes
lon_spacing = 0:30:330;
[txlats, txlons] = repeat_crd(txlats, txlons, lon_spacing);
[rxlats, rxlons] = repeat_crd(rxlats, rxlons, lon_spacing);


%% Load and reformat SAMI

for t = 1:length(times)
    time = times(t);
    sami = loadstruct(filename(sami_fn_fmt, time));

    lons = sami.lon;
    lons(lons > 180) = lons(lons > 180) < 360;

    %% Loop through tx/rx locations
    for l1 = 1:length(txlats)
        %% Generate oblique ionograms and save
        txloc = [txlats(l1), txlons(l1), alt];
        rxloc = [rxlats(l1), rxlons(l1), alt];

        if txlons(l1) > max(lons) || rxlons(l1) > max(lons) || ...
                txlons(l1) < min(lons) || rxlons(l1) < min (lons)
            fprintf('Transmitter or receiver longitude outside grid - skipping\n')
            continue
        end
        homed_rays = gen_ionogram(freqs, OX_mode, txloc, rxloc, ...
            sami.iono_en_grid, sami.iono_en_grid, sami.collision_freq, ...
            sami.iono_grid_parms, sami.geomag_grid_parms, sami.Bx, sami.By, sami.Bz, ...
            maxdist, tol, homing_tol_m);

        [frq, rg] = calc_ionogram(homed_rays);

        if isempty(frq)
            fprintf('No valid rays - skipping\n')
            continue
        end
        % save
        out_fn = sprintf(filename(out_fn_fmt, time), ...
            txloc(1), txloc(2), rxloc(1), rxloc(2), txloc(3));
        savestruct(out_fn, homed_rays)
        fprintf('Saved to %s\n', out_fn)

    end

end



%% coordinate meshing function
function [lats, lons] = repeat_crd(lats, lons, lon_spacing)
loninc = repmat(lon_spacing, [length(lats), 1]);
lats = repmat(lats, [1, length(loninc)]);
lons = repmat(lons, [1, length(loninc)]);
lons = lons + loninc;
lats = lats(:);
lons = lons(:);
lons(lons >= 180) = lons(lons >= 180) - 360;

return
end



