function [Bx, By, Bz] = gen_bfield(lats, lons, hts, yr)
%% generate a geomagnetic field for PHaRLAP from IGRF
% % e.g. 
% lats = 40:5:90; 
% lons = 0:2; 
% hts = 0:100:600; 
% yr = 2012;
%
% cartv = sphcart([64E5+ht_3d_m(:), deg2rad(lat_3d(:)), deg2rad(lon_3d(:))]);
% X = reshape(cartv(:, 1), size(Bx));
% Y = reshape(cartv(:, 2), size(Bx));
% Z = reshape(cartv(:, 3), size(Bx));
% 
% plot3(squeeze(X(:, 1, :)), squeeze(Y(:, 1, :)), squeeze(Z(:, 1, :)), '.')
% 
% quiver3(squeeze(X(:, 1, :)), squeeze(Y(:, 1, :)), squeeze(Z(:, 1, :)),...
%     squeeze(Bx(:, 1, :)), squeeze(By(:, 1, :)), squeeze(Bz(:, 1, :)), ...
%     'ShowArrowHead','off')
% 
% contourf(lats, hts, squeeze(B(:, 1, :))')

%%
[lat_3d, lon_3d, ht_3d_m] = ndgrid(lats, lons, hts * 1E3);
NEU = igrfmagm(ht_3d_m(:), lat_3d(:), lon_3d(:), ones(size(ht_3d_m(:))) .* yr) * 1E-9;
NEU = reshape(NEU,[size(lat_3d,1),size(lat_3d,2),size(lat_3d,3),3]);
[Bx, By, Bz] = enu2ecefv(...
    NEU(:, :, :, 2), NEU(:, :, :, 1), -NEU(:, :, :, 3), lat_3d, lon_3d);

B = sqrt(Bx.^2 + By.^2 + Bz.^2);




