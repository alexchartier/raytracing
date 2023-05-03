
%%  inputs
fn_fmt = 'data/wspr/links/{YYYY-mmm-dd-HHMM}.nc';
latbins = -90:10:90;
ltbins = 0:24;
times = datenum(2019, 3, 2):1/24:datenum(2019, 4, 1);

%% load
% data_full = read_netcdf(fn(2).name);  % skip day 1 (spinup)
clear data_full
for t = 1:length(times)
    fn = filename(fn_fmt, times(t));

    % load
    try
        data = read_netcdf(fn);
    catch
        continue
    end

    data.time = ones(size(data.home)) .* times(t);
    fieldn = fieldnames(data);

    if t == 1
        data_full = data;
    else
        for f = fieldn'
            data_full.(f{1}) = [data_full.(f{1}); data.(f{1})];
        end
    end
end


%% clean up/derive relevant params
npts = size(data_full.home);
data_full.lat = ones(npts) .* NaN;
data_full.lon = ones(npts) .* NaN;
for l = 1:npts
    % Great circle calcs
    [ilat, ilon] = ...
        interpm([data_full.txlocs(l, 1), data_full.rxlocs(l, 1)], ...
        [data_full.txlocs(l, 2), data_full.rxlocs(l, 2)], 360, 'gc');
    data_full.lat(l) = median(ilat);
    data_full.lon(l) = median(ilon);
end

dt = datenum(data_full.time);
lt = (data_full.lon / 360 + (dt- floor(dt))) * 24;
lt(lt < 0) = lt(lt < 0) + 24;
lt(lt >= 24) = lt(lt >= 24) - 24;
data_full.lt = lt;



%% binned histograms
ll = mean([latbins(1:end-1); latbins(2:end)]);
laterrs = zeros(length(ll), 4);
for l = 1:length(latbins) - 1
    lati = data_full.lat > latbins(l) & data_full.lat <= latbins(l+1);
    laterrs(l, 1) = sum(lati);
    laterrs(l, 2) = sum(data_full.home(lati));
end

hh = mean([ltbins(1:end - 1); ltbins(2:end)]);
lterrs = zeros(length(hh), 4);
for l = 1:length(ltbins) - 1
    lti = data_full.lt > ltbins(l) & data_full.lt <= ltbins(l+1);
    lterrs(l, 1) = sum(lti);
    lterrs(l, 2) = sum(data_full.home(lti));
end


%% day/night comparison
dayi = data_full.lt > 6 & data_full.lt < 18;
nighti = ~dayi;
daypct = 100 * sum(data_full.home(dayi)) ./ sum(dayi);
nightpct = 100 * sum(data_full.home(nighti)) ./ sum(nighti);


%% plot
fs = 22;
subplot(2, 1, 1)

bar(ll, laterrs)
xlabel('Midpoint Latitude (deg)')
ylabel('# WSPR links')
legend({"Reported", "Predicted by SAMI3"})
text(-80, 4E4, sprintf('%i total links\n%i predicted\n%1.1f%% success', ...
    length(data_full.home), sum(data_full.home), ...
    100 * sum(data_full.home)./ length(data_full.home)), "FontSize", fs ...
    )

set(gca, 'FontSize', fs)

grid on
grid minor

subplot(2, 1, 2)
bar(hh, lterrs)
xlabel('Midpoint Local Time (hour)')
ylabel('# WSPR links')

legend({"Reported", "Predicted by SAMI3"})
text(1, 7E3, sprintf('6 - 18 LT: %1.1f%% success \n18 - 6 LT: %1.1f%% success\n', ...
    daypct, nightpct), "FontSize", fs ...
    )

set(gca, 'FontSize', fs, 'ylim', [0, 8000])

grid on
grid minor



