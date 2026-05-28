%% prep_iri_for_raytracing.m

%% calculate IRI equivalent to SAMI3 structure
% TODO: FIXME (order is wrong)
sami_fn_fmt = '/Users/chartat1/data/sami3/2017_tid/sami_mat/{YYYY-mm-dd_HHMM}.mat';
iri_fn_fmt = '/Users/chartat1/data/iri/{YYYY-mm-dd_HHMM}.mat';
time = datenum(2017, 1, 10, 18, 0, 0);% :5/60/24:datenum(2017, 1, 11, 0, 0, 0);
sami = loadstruct(filename(sami_fn_fmt, time));

alti = sami.alt <= alt_cutoff;
sami.alt = sami.alt(alti);
sami.dene = sami.dene(alti, :, :);

R12 = 28;
lon = sami.lon;
lon(lon>180) = lon(lon>180) - 360;
lon = sort(lon);
[iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    Bx, By, Bz, geomag_grid_parms] = gen_iono_geomag_grids(...
    sami.alt, sami.lat, sami.lon, time, R12, sami.geomag_grid_parms(8), ...
    sami.geomag_grid_parms(2), sami.geomag_grid_parms(5));

iri.time = time;
iri.alt = sami.alt;
iri.lat = sami.lat;
iri.lon = sami.lon;
dene = permute(iono_en_grid, [3, 1, 2]);
dene_reorder = cat(3, dene(:, :, iri.lon > 180), dene(:, :, iri.lon <= 180));
iri.dene = dene_reorder; 
iri.iono_en_grid = iono_en_grid;
iri.collision_freq = collision_freq;
iri.Bx = Bx; 
iri.By = By; 
iri.Bz = Bz; 
iri.iono_grid_parms = sami.iono_grid_parms;
iri.geomag_grid_parms = sami.geomag_grid_parms;

savestruct(filename(iri_fn_fmt, time), iri)
