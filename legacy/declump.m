function [latI, lonI, idx] = declump(lat, lon, npts_wanted)
%% declump.m
% Function to de-clump data, downsampling to a required number of points
% NOTE: we are not doing this perfectly or accounting for the longitude
% wrap
%
% % Example
% 
% lats = [-90:90];
% lons = [-180:2:180];
% rp1 = randperm(length(lats));
% rp2 = randperm(length(lons));
% npts = 50;
% lats = lats(rp1(1:npts));
% lons = lons(rp2(1:npts));
% npts_wanted = 20;
% [latI, lonI, idx] = declump(lats, lons, npts_wanted);
% hold on
% scatter(lons, lats, 'b')
% scatter(lonI, latI, '.k')
% xlim([-180, 180])
% ylim([-90, 90])
% hold off



%% Calculate distances between them
x = [lat; lon]';
[~, dists] = knnsearch(x, x, 'K', 10);
distarr = sum(dists, 2);

[B, idx1] = maxk(distarr, floor(npts_wanted * 0.7));

nvals = 1:size(x, 1);
idx2 = nvals(~ismember(nvals, idx1));
rp = randperm(length(idx2));
idx = [idx1; idx2(rp(1:ceil(npts_wanted * 0.3)))'];


%% Select the n links you want
latI = lat(idx);
lonI = lon(idx);
