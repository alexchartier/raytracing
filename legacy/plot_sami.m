%% Set filenames
fname = '/Users/chartat1/data/sami3/2019/sami3_regulargrid_elec_density_2019Mar13.nc';
ncdisp(fname)


%% load
dene = ncread(fname, 'dene0');
lats =  ncread(fname, 'lat');
lons =  ncread(fname, 'lon');
alts =  ncread(fname, 'alt');
times =  ncread(fname, 'time');


%% plot
t = 1;
dene_t = dene(:, :, :, t);
clf
x0=10;
y0=10;
width=1600;
height=1200;

set(gcf,'position',[x0,y0,width,height])
tiledlayout(2, 1, 'TileSpacing','compact');

nexttile
contourf(lats, alts, squeeze(dene_t(:, 1, :)))
ylim([0, 600])
ylabel('Alt (km)')
xlabel('Lat (deg)')
title(fname)
set(gca, 'FontSize', 24)


nexttile
alti = alts == max(alts(alts <350));
contourf(lons, lats, squeeze(dene_t(alti, :, :))')
xlabel('Lon (deg)')
ylabel('Lat (deg)')
set(gca, 'FontSize', 24)








































