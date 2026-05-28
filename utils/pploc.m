function pplocs = pploc(ray, rxloc)
%% calculate the lat/lon at the piercepoint(s) of the alt shell
nrays = length(ray);
P_lat = ones(nrays, 1) * nan;
P_lon = ones(nrays, 1) * nan;
for i = 1:nrays
    upidx = diff(ray(i).height) > 0;
    if sum(upidx) <= 1  % ray went through to the ground
        P_lat(i) = NaN;
        P_lon(i) = NaN;
    else
        P_lat(i) = interp1(ray(i).height(upidx), ray(i).lat(upidx), rxloc(3));
        P_lon(i) = interp1(ray(i).height(upidx), ray(i).lon(upidx), rxloc(3));
    end
end
pplocs = [P_lat, P_lon];
