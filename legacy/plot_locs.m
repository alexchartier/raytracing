%% plot the data locations for the paper

D = loadstruct('wspr_locs.mat');

wspr_t = D.wspr_t;
wspt_out = D.wspr_out;

lons_iono = D.lons;
lons_iono(lons_iono > 180)  = lons_iono(lons_iono > 180) - 360;
lats_iono = D.lats;


slist_full = asciiread('data/gps/sitelist.txt');
slist_short = asciiread('data/gps/sitelist_150.txt');
latlon_full = str2num(slist_full(:, 8:end));
latlon_short = str2num(slist_short(:, 8:end));

%%
clf
hold on
m_proj('Miller Cylindrical', 'latitude', [-90 90], 'longitude', [-180, 180])

m_plot(latlon_full(:, 2), latlon_full(:, 1), 'og', 'MarkerSize', 10)
m_plot(latlon_short(:, 2), latlon_short(:, 1), '.g', 'MarkerSize', 40)


m_plot(wspr_t.lonI, wspr_t.latI, 'ob', 'MarkerSize', 10)
m_plot(wspr_out.lonI, wspr_out.latI, '.b', 'MarkerSize', 40)


m_plot(lons_iono, lats_iono, '.r', 'MarkerSize', 40)

m_grid('color', 'k', 'FontSize', 20)
m_coast('color', 'k');

hold off
legend({'', 'GPS stations (all)', 'GPS stations (selected)', ...
    'WSPR midpoints (all)', 'WSPR midpoints (selected)', 'Ionosonde Stations'}, 'FontSize', 24)
