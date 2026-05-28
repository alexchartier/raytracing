%% inputs
UT = [2012, 1, 1, 0, 0];
R12 = 100;
model.lat = -80:5:80;
model.lon = 0:15:345;
model.alt = 100:5:2095;
mod_fn_fmt = '/Users/chartat1/data/nebula/model_input/iri_2000km/{YYYY-mm-dd_HHMM}.mat';

%% generate foEs map;
model.time = datenum(UT(1), UT(2), UT(3), UT(4), UT(5), 0);
ht_min = min(model.alt);
ht_inc = model.alt(2) - model.alt(1);
num_heights = length(model.alt);

% foEs = (sin(-2:length(model.lat) - 3) + 2) * 3;

% foEs = repmat(foEs, [length(model.lon), 1])';
% hmE = foEs + 100;
model.dene = zeros([length(model.alt), length(model.lat), length(model.lon)]);

% call IRI2020
for l1 = 1:length(model.lat)
    for l2 = 1:length(model.lon)
        % N = calc_chapman_prof(model.alt, freq2elec(foEs(l1, l2) * 1E6), ...
        %     hmE(l1, l2), 10 - foEs(l1,l2)/2, 1);

        [iono, iono_extra] = iri2016(model.lat(l1), model.lon(l2), R12, UT, ...
            ht_min, ht_inc, num_heights);
        prof = iono(1, :);%   + N;
        model.dene(:, l1, l2) = prof / 1E6;
    end
end

%% Generage geomagnetic field
[iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    Bx, By, Bz, geomag_grid_parms] = gen_grid_parms(model);
model.iono_en_grid = iono_en_grid;
model.collision_freq = collision_freq;
model.iono_grid_parms = iono_grid_parms;
model.geomag_grid_parms = geomag_grid_parms;
model.Bx = Bx;
model.By = By;
model.Bz = Bz;

%% Save
savestruct(filename(mod_fn_fmt, model.time), model)
fprintf('Saved to %s\n', filename(mod_fn_fmt, model.time))