function D = load_sami(sami_fn, vars, stime)
%% Function load_sami
% sami_fn = 'test.nc';
% 
% vars = {'alt', 'lat', 'lon', 'dene0', 'time'};
% D = load_sami(sami_fn, vars);
% contourf(D.lon, D.lat, squeeze(D.dene0(30, :, :, 100))')

%%

D = [];
for v = vars
    D.(v{1}) = double(ncread(sami_fn, v{1}));
end

D.time = stime + D.time / 24;


