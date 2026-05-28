function [latMid, lonMid] = midpointLatLon(lat1, lon1, lat2, lon2)
% midpoint of two lat long cord on a sphere, all units are deg
Bx = cosd(lat2) * cosd(lon2-lon1);
By = cosd(lat2) * sind(lon2-lon1);
latMid = atan2d(sind(lat1) + sind(lat2), ...
               sqrt( (cosd(lat1)+Bx)*(cosd(lat1)+Bx) + By*By ) );
lonMid = lon1 + atan2d(By, cosd(lat1) + Bx);