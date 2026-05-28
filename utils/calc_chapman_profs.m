function iono_en_grid = calc_chapman_profs(lats, lons, alts, Np, hp, H, c)
%% Create a Chapman ionosphere over a 3D grid
% h - altitudes
% Np - peak density
% hp - peak height
% H - scale height
% c - chapman parameter (typically 0 - 1
%
% lats = 40:2:60;
% lons = -75:5:-65;
% alts = 150:2:600;
% Np = length(lats).* 5E5;
% hp = length(lats).* 300;
% H = length(lats).* 100;
% c = 1;
%
% iono_en_grid = calc_chapman_profs(lats, lons, alts, Np, hp, H, c);
% plot(squeeze(iono_en_grid(1, 1, :)), alts)


%%
iono_en_grid = zeros([length(lats), length(lons), length(alts)]);
for i = 1:length(lats)
    for j = 1:length(lons)
        iono_en_grid(i, j, :) = ...
            calc_chapman_prof(alts, Np(i, j), hp(i, j), H(i, j), c);
    end
end
