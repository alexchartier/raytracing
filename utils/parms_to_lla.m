function [lat, lon, ht] = parms_to_lla(parms)
%% parms_to_lla
% convert the PHaRLAP parameters to lat/lon/alt

%%
lat = arange(parms(1), parms(2), parms(3));
lon = arange(parms(4), parms(5), parms(6));
ht = arange(parms(7), parms(8), parms(9));
