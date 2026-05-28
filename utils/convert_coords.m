function loc = convert_coords(XYZ)
%% calculate the lat(°)/lon(°)/alt(km) location from ECEF XYZ (m)

RE = 6371E3;

loc_RLL_rad = cartsph(XYZ);
loc = [rad2deg(loc_RLL_rad(2)), rad2deg(loc_RLL_rad(3)), (loc_RLL_rad(1) - RE) / 1E3];

if loc(2) > 180
    loc(2) = loc(2) - 360;
end

