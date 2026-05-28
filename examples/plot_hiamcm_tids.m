%% plot_hiamcm_tids.m
% Plot the TIDs from the SAMI3/HIAMCM run

%%

in_fn_fmt = '/Users/chartat1/data/sami3/2017_tid/sami_mat/{YYYY-mm-dd_HHMM}.mat';
plot_dirn = '/Users/chartat1/data/sami3/2017_tid/plots/';
alt_time_plt = [plot_dirn, 'alt_time_{YYYY-mm-dd}.png'];
latlon_mov_plt = [plot_dirn, 'mov/latlon_{YYYY-mm-dd-HHMM}.png'];
altlat_mov_plt = [plot_dirn, 'mov_local/altlat_{YYYY-mm-dd-HHMM}.png'];

alt_lat_plt = [plot_dirn, 'alt_lat_{YYYY-mm-dd-HHMM}.png'];
times = datenum(2017, 1, 10, 0, 0, 0):5/60/24:datenum(2017, 1, 10, 23, 55, 0);

lon = 282.5;
lat = 40;

%% load
D = loadstruct(filename(in_fn_fmt, times(1)));
ne_grid = zeros(length(D.alt), length(times));
for t = 1:length(times)
    D = loadstruct(filename(in_fn_fmt, times(t)));
    ne_grid(:, t) = D.dene(:, D.lat == lat, round(D.lon*10) == round(lon*10));
end

%% plot alt vs time

contourf(times, D.alt, ne_grid);
title(sprintf('SAMI3/HIAMCM Electron Density, %1.1f N, %1.1f E', lat, lon))
xlabel('Time (UT)')
ylabel('Alt (km)')
colorbar
datetick('keeplimits')
export_fig(filename(alt_time_plt, times(1)))

%% Plot lat/lon contours at multiple times
figure('Position', [800, 600, 800, 600])
idx = 105;
for t = 1:length(times)
    %%
    D = loadstruct(filename(in_fn_fmt, times(t)));
    [~, hC] = contourf(D.lon, D.lat, elec2freq(squeeze(D.dene(idx, :, :))) / 1E3, 100);
    set(hC, 'LineStyle', 'none')
    clim([0, 12])
    a =colorbar;
    a.Label.String = 'Plasma Freq (MHz)';
    xlabel('Lon (°)')
    ylabel('Lat (°)')
    title(sprintf('%s, %i km alt.', filename('{YYYY-mm-dd HH:MM UT}', times(t)), D.alt(idx)))
    pause(0.1)
    export_fig(filename(latlon_mov_plt, times(t)))

end


%% Plot lat vs alt
times = datenum(2017, 1, 10, 12, 0, 0):5/60/24:datenum(2017, 1, 10, 23, 55, 0);
lats = lat - 10:2:lat+10;
lons = 282.5;
close
figure('Position', [800, 600, 800, 600])


for t  = 1:length(times)
    D = loadstruct(filename(in_fn_fmt, times(t)));
    dene = squeeze(D.dene(:, ismember(D.lat, lats), ismember(D.lon, lons)));
    [~, hC] = contourf(lats, D.alt, elec2freq(dene) / 1E3);
    xlabel('Lat (°)')
    ylabel('Alt (km)')
    set(hC, 'LineStyle', 'none')
    clim([0, 6])
    a =colorbar;
    a.Label.String = 'Plasma Freq (MHz)';
    title(sprintf('%s, %1.1f° E', filename('{YYYY-mm-dd HH:MM UT}', times(t)), lons))
    pause(0.1)
    export_fig(filename(altlat_mov_plt, times(t)))
end






















