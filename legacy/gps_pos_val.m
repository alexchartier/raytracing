%% gps_pos_val.m
% plot the box-whiskers

%% filenames
fn = rdir('~/data/gps_pos_errs/zero_init_pos_errs*.nc');
latbins = -90:10:90;
ltbins = 0:24;


%% load
data_full = read_netcdf(fn(2).name);  % skip day 1 (spinup)
fieldn = fieldnames(data_full);
for i = 3:length(fn)
    % load

    data = read_netcdf(fn(i).name);
    for f = fieldn'
        data_full.(f{1}) = [data_full.(f{1}); data.(f{1})];
    end
end

good_idx = data_full.corr_pos_err < 50;
for f = fieldn'
    data_full.(f{1}) = data_full.(f{1})(good_idx);
end

data_full.time = datetime(data_full.time, 'ConvertFrom', 'posixtime');
dt = datenum(data_full.time);

lt = (data_full.lon / 360 + (dt- floor(dt))) * 24;
lt(lt < 0) = lt(lt < 0) + 24;
lt(lt >= 24) = lt(lt >= 24) - 24;
data_full.lt = lt;


% %% process site-by-site
% 
% % set up the holder
% sites = unique(data_full.site);
% for s = 1:length(sites)
%     sites{s} = sprintf('a%s', sites{s});
% end
% 
% sites = cellstr(sites');
% f = sites;
% f{2, 1} = [];
% D = struct(f{:});
% 
% 
% % store everything
% for s = sites
%     sn = s{1};
%     si = data_full.site == sn(2:end);
%     D.(sn).lat = unique(data_full.lat(si));
%     D.(sn).lon = unique(data_full.lon(si));
%     D.(sn).raw_pos_err = data_full.raw_pos_err(si);
%     D.(sn).corr_pos_err = data_full.raw_pos_err(si);
%     D.(sn).lt = data_full.lt(si);
% end


%% binned histograms
ll = mean([latbins(1:end-1); latbins(2:end)]);
laterrs = zeros(length(ll), 4);
for l = 1:length(latbins) - 1
    lati = data_full.lat > latbins(l) & data_full.lat <= latbins(l+1);
    laterrs(l, 1) = median(data_full.raw_pos_err(lati));
    laterrs(l, 2) = max(data_full.raw_pos_err(lati));
    laterrs(l, 3) = median(data_full.corr_pos_err(lati));
    laterrs(l, 4) = max(data_full.corr_pos_err(lati));
end

hh = mean([ltbins(1:end - 1); ltbins(2:end)]);
lterrs = zeros(length(hh), 4);
for l = 1:length(ltbins) - 1
    lti = data_full.lt > ltbins(l) & data_full.lt <= ltbins(l+1);
    lterrs(l, 1) = median(data_full.raw_pos_err(lti));
    lterrs(l, 2) = max(data_full.raw_pos_err(lti));
    lterrs(l, 3) = median(data_full.corr_pos_err(lti));
    lterrs(l, 4) = max(data_full.corr_pos_err(lti));
end


%% day/night comparison
dayi = data_full.lt > 6 & data_full.lt < 18;
nighti = ~dayi;
highlati = abs(data_full.lat) > 60;
midlati = abs(data_full.lat) < 60| abs(data_full.lat) > 30;
lowlati = abs(data_full.lat) < 30;

errpct = 100 * (1 - data_full.corr_pos_err ./ data_full.raw_pos_err);

daypct = median(errpct(dayi));
nightpct = median(errpct(nighti));
    
highpct = median(errpct(highlati));
midpct = median(errpct(midlati));
lowpct = median(errpct(lowlati));


%% plot
clf
fs = 22;

subplot(2, 1, 1)

bar(ll, laterrs(:, 1:2:3))
xlabel('Latitude (deg)')
ylabel('Median 3D GPS position error (m)')
legend({"Uncorrected", "Corrected by SAMI3"})
text(-80, 4.5, sprintf('Overall median & max raw err: %1.1f & %1.1f m\nOverall median & max corr. err: %1.1f & %1.1f m',...
    median(data_full.raw_pos_err), max(data_full.raw_pos_err),  ...
    median(data_full.corr_pos_err), max(data_full.corr_pos_err)), ...
    "FontSize", fs)
text(-80, 3, sprintf('> 60 deg: %1.1f%% improvement \n30 - 60 deg: %1.1f%% improvement\n< 30 deg: %1.1f%% improvement\n', ...
    highpct, midpct, lowpct), "FontSize", fs ...
    )

set(gca, 'FontSize', fs)

grid on
grid minor

subplot(2, 1, 2)
bar(hh, lterrs(:, 1:2:3))
xlabel('Local Time (hour)')
ylabel('Median 3D GPS position error (m)')
legend({"Uncorrected", "Corrected by SAMI3"})
text(1, 4, sprintf('6 - 18 LT: %1.1f%% improvement \n18 - 6 LT: %1.1f%% improvement\n', ...
    daypct, nightpct), "FontSize", fs ...
    )

set(gca, 'FontSize', fs)

grid on
grid minor


% %% box/whisker plots
% 
% hold on
% 
% for s = sites
%     sn = s{1};
%     % plot the box/whiskers
%     lat = D.(sn).lat(1);
%     errs = D.(sn).corr_pos_err;
%     boxplot(errs, 'positions', lat, 'labels', lat, 'symbol', '')
% end
% ylim([0, 15])
% ylabel('3D GPS position error (m)')
% xlabel('latitude')
% hold off




