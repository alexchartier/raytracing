function n = interp_3d(lats, lons, alts, n_3d, loc)
%% 3D stepwise interpolation of density field in lla order

lati = loc(1); 
loni = loc(2); 
alti = loc(3);
a1 = max(alts(alts <= alti));
a2 = min(alts(alts >= alti));
n1 = interp2(lats, lons, squeeze(n_3d(:, :, alts == a1))', lati, loni);
n2 = interp2(lats, lons, squeeze(n_3d(:, :, alts == a2))', lati, loni);

if n1 == n2
    n = n1;
else
    n = interp1([a1, a2], [n1, n2], alti);
end