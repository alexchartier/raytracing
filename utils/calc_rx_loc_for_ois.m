function rxloc = calc_rx_loc_for_ois(txloc, refloc, alongtrack_spacing_km)
%% 
% rxloc = calc_rx_loc_for_ois(txloc, refloc, alongtrack_spacing_km)
% 
% reckons a location the specified distance away from txloc in the
% direction of refloc, at the same altitude as both
wgs84 = wgs84Ellipsoid;
RE = mean([wgs84.SemimajorAxis, wgs84.SemiminorAxis]);
arclen = 360 * alongtrack_spacing_km * 1E3 / (2 * pi * (RE + txloc(3) * 1E3));
az = azimuth(txloc(1), txloc(2), refloc(1), refloc(2), wgs84Ellipsoid, 'degrees');
[rxlat, rxlon] = reckon(txloc(1), txloc(2), arclen, az);
rxloc = [rxlat, rxlon, txloc(3)];


