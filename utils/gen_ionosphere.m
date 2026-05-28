function [iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    Bx, By, Bz, geomag_grid_parms, lats, lons, alts] = ...
    gen_ionosphere(sat, time, timeidx, R12, out_iono, B_ht_inc, B_lat_inc, B_lon_inc, ...
    I_ht_inc, I_lat_inc, I_lon_inc, pad_lat, pad_alt, load_iono)

% Load/compute ionosphere
if load_iono
    struc = load(join([out_iono, 'iono_en_grid.mat']), '-mat', 'iono_en_grid');
    iono_en_grid = struc.iono_en_grid;
    iono_en_grid_5 = iono_en_grid;
    struc = load(join([out_iono, 'collision_freq.mat']));
    collision_freq = struc.collision_freq;
    struc = load(join([out_iono, 'iono_grid_parms.mat']));
    iono_grid_parms = struc.iono_grid_parms;
    struc = load(join([out_iono, 'Bx.mat']));
    Bx = struc.Bx;
    struc = load(join([out_iono, 'By.mat']));
    By = struc.By;
    struc = load(join([out_iono, 'Bz.mat']));
    Bz = struc.Bz;
    struc = load(join([out_iono, 'geomag_grid_parms.mat']));
    geomag_grid_parms = struc.geomag_grid_parms;
    struc = load(join([out_iono, 'lons']));
    lons = struc.lons;
    struc = load(join([out_iono, 'lats']));
    lats = struc.lats;
    struc = load(join([out_iono, 'alts']));
    alts = struc.alts;
else
    % find maximum alt,lat,lon
    max_alt = max(sat.alt{:,:},[],"all");
    min_lat = min(sat.lat{:,:},[],"all");
    max_lat = max(sat.lat{:,:},[],"all");
    
    lons = -180:I_lon_inc:180;
    alts = 92:I_ht_inc:max_alt + pad_alt;
    if min(min_lat) - pad_lat >= -90
        start_lat = min_lat - pad_lat;
    else
        start_lat = -90;
    end
    if max_lat + pad_lat <= 90
        end_lat = max_lat + pad_lat;
    else
        end_lat = 90;
    end
    lats = start_lat:I_lat_inc:end_lat;

    [iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
        Bx, By, Bz, geomag_grid_parms] = gen_iono_geomag_grids( ...
        alts, lats, lons, time, R12, B_ht_inc, B_lat_inc, B_lon_inc);
    save(join([out_iono, 'iono_en_grid']), 'iono_en_grid');
    save(join([out_iono, 'iono_en_grid_5']), 'iono_en_grid_5');
    save(join([out_iono, 'collision_freq']), 'collision_freq');
    save(join([out_iono, 'iono_grid_parms']), 'iono_grid_parms');
    save(join([out_iono, 'Bx']), 'Bx');
    save(join([out_iono, 'By']), 'By');
    save(join([out_iono, 'Bz']), 'Bz');
    save(join([out_iono, 'geomag_grid_parms']), 'geomag_grid_parms');
    save(join([out_iono, 'lons']), 'lons');
    save(join([out_iono, 'lats']), 'lats');
    save(join([out_iono, 'alts']), 'alts');
end
