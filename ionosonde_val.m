%% ionosonde_val.m
% plot the box-whiskers

%% filenames
fn = rdir('data/ionosonde_val/iono*.nc');

keys = {'nmF2', 'nmF2_sami', 'lat', 'lt'};


%%
clf
hold on

for i = 1:length(fn)
    % load
    for k = keys
        data.(k{1}) = ncread(fn(i).name, k{1});
        
    end
    if isempty(data.nmF2)
        continue
    end

    % plot
    lat = data.lat(1);
        errs = data.nmF2_sami - data.nmF2;
        boxplot(errs, 'positions', lat, 'labels', lat, 'symbol', '')
end

xlabel('latitude')
hold off


%     nmf2 = ncread(fn(i).name, 'nmF2');
%     nmf2_sami = ncread(fn(i).name, 'nmF2_sami');
%     lt = ncread(fn(i).name, 'lt');