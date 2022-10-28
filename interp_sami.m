function D2 = interp_sami(D, alts, time)
% 
% D = load_sami('test.nc');
% alts = 90:2:400;
% time = datetime(2014, 5, 23, 6, 0, 0);
% D2 = interp_sami(D, alts, time)
%  
% % plot
% contourf(sami_t.lon, sami_t.lat, sami_t.dene(:, :, 100))

%%
D2 = D;
size_D = size(D.dene0);
size_D(1) = length(alts);
D2.dene = nan(size_D(1:end - 1));
D2.alt = alts;
D2 = rmfield(D2, 'dene0');

tind = D.time == time;

dene_in = D.dene0(:, :, :, tind);

for i1 = 1:length(D.lon)
    for i2 = 1:length(D.lat)
        D2.dene(:, i1, i2) = interp1(D.alt, dene_in(:, i1, i2), alts);
    end
end
D2.dene = permute(D2.dene, [3, 2, 1]);

%% Correct to -180/180 from 0/360 and wrap in lat and lon

if max(D2.lon) >= 180
    lonidx = D2.lon >= 180;
    D2.lon = [D2.lon(lonidx) - 360; D2.lon(~lonidx); 180];
    D2.dene = cat(2, D2.dene(:, lonidx, :), D2.dene(:, ~lonidx, :));
    D2.dene = cat(2, D2.dene, D2.dene(:, 1, :));
end



































