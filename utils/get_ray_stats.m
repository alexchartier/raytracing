function [maxne, minht, nmf2, hmf2, loc] = get_ray_stats(rays, iono_en_grid, lats, lons, alts)


% Check that we have a 'good' golden ray
minht = min(rays{1}.height);
loci = rays{1}.height == minht;
maxne = rays{1}.electron_density(loci);
loc = [rays{1}.lat(loci), rays{1}.lon(loci), rays{1}.height(loci)];
nmf2_grid = max(iono_en_grid, [], 3);
hmf2_grid = zeros(size(nmf2_grid));

nmf2 = interp2(lats,lons, nmf2_grid', loc(1), loc(2));
for i = 1:length(lats)
    for j = 1:length(lons)
        hmf2_grid(i, j) = alts(iono_en_grid(i, j, :) == nmf2_grid(i,j));
    end
end

hmf2 = interp2(lats,lons, hmf2_grid', loc(1), loc(2));

fprintf('fit / input / quantity\n%1.1e %1.1e NmF2\n%1.0f %1.0f hmf2\n', maxne, nmf2, minht, hmf2)
