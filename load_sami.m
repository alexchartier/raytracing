function D = load_sami(sami_fn)
%% Function load_sami
% sami_fn = 'test.nc';
% D = load_sami(sami_fn);
% contourf(D.lon, D.lat, squeeze(D.dene0(30, :, :, 100))')

%%
vars = {'alt', 'lat', 'lon', 'dene0', 'time'};
D = [];
for v = vars
    D.(v{1}) = double(ncread(sami_fn, v{1}));
end

D.time = datetime(D.time, 'ConvertFrom', 'posixtime');


