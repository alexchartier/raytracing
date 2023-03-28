%% ionosonde_val.m
% plot the box-whiskers

%% filenames
fn = rdir('data/ionosonde_val/iono*.nc');
stime = datenum(2019, 3, 2);
latbins = -90:10:90;
ltbins = 0:24;


%% Load
data_full = read_netcdf(fn(1).name);  % skip day 1 (spinup)

fieldn = fieldnames(data_full);
for i = 2:length(fn)
    % load

    data = read_netcdf(fn(i).name);
    for f = fieldn'
        data_full.(f{1}) = [data_full.(f{1}); data.(f{1})];
    end
end

data_full.time = datetime(data_full.POSIXtime, 'ConvertFrom', 'posixtime');
dt = datenum(data_full.time);
tidx = dt >= stime;
finidx = ~isnan(data_full.nmF2_sami);
idx = tidx & finidx;
for f = fieldn'
    data_full.(f{1}) = data_full.(f{1})(idx);
end
data_full.time = data_full.time(idx);
dt = dt(idx);

lt = (data_full.lon / 360 + (dt- floor(dt))) * 24;
lt(lt < 0) = lt(lt < 0) + 24;
lt(lt >= 24) = lt(lt >= 24) - 24;
data_full.lt = lt;

data_full.errs = data_full.nmF2_sami - data_full.nmF2;

%% box-whisker plots
clf
% Lat
subplot(2, 1, 1)
hold on
ll = mean([latbins(1:end-1); latbins(2:end)]);
laterrs = zeros(length(ll), 4);
for l = 1:length(latbins) - 1
    lati = data_full.lat > latbins(l) & data_full.lat <= latbins(l+1);
    d = data_full.errs(lati);
    b = boxchart(ll(l) * ones(size(d)), d, 'BoxWidth', 10, 'MarkerStyle', '.');
    b.JitterOutliers = 'on';

end
ylabel('SAMI3 NmF2 err. (el. m-3)')
xlabel('Latitude (10 deg bins)')
xlim([-90, 90])
plot(latbins, zeros(size(latbins)), '-k')
set(gca, 'FontSize', 24)
grid on
grid minor

text(50, -9E11, sprintf(...
    'Overall NmF2 errors\nMedian: %1.1f x10^{11} el. m^{-3}\nMin: %1.1f x10^{11} el. m^{-3}\nMax: %1.1f x10^{11} el. m^{-3}\n', ...
    median(data_full.errs)/1E11, min(data_full.errs)/1E11, max(data_full.errs)/1E11), "FontSize", fs)

% LT
subplot(2, 1, 2)
hold on
hh = mean([ltbins(1:end - 1); ltbins(2:end)]);
lterrs = zeros(length(hh), 4);
for l = 1:length(ltbins) - 1
    lti = data_full.lt > ltbins(l) & data_full.lt <= ltbins(l+1);
    d = data_full.errs(lti);
    b = boxchart(hh(l) * ones(size(d)), d, 'BoxWidth', 1, 'MarkerStyle', '.');
    b.JitterOutliers = 'on';
end
ylabel('SAMI3 NmF2 err. (el. m-3)')
xlabel('Local Time (1 hr bins)')
xlim([0, 24])
plot(ltbins, zeros(size(ltbins)), '-k')
set(gca, 'FontSize', 24)
grid on
grid minor






% %% siteplot
% clf
% m_proj('Miller Cylindrical', 'latitude', [-80 80], 'longitude', [0, 360])
% 
% m_plot(lons, lats, '.r', 'MarkerSize', 40)
% 
% 
%     m_grid('color', 'k', 'FontSize', 20)
%     m_coast('color', 'k');
% 
%     legend({'', 'Ionosonde Stations'}, 'FontSize', 24)
% 
% 
% %     nmf2 = ncread(fn(i).name, 'nmF2');
% %     nmf2_sami = ncread(fn(i).name, 'nmF2_sami');
% %     lt = ncread(fn(i).name, 'lt');






