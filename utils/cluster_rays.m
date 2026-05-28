function [minima, close_ray_mat] = cluster_rays(elvarr, azarr, rays, rxloc, maxdist)

%% cluster the rays
close_rays = [];
ids = [];
min_di = 1E9;

% TODO First go through and kick out the 'bad' ones'

% for r = 1:length(rays_in)
%     disp('r')
% end

%% 

for r = 1:length(rays)
    [di, close_pt, id] = ray_dist(rays(r), rxloc);
    % disp(di)
    if di < maxdist
        rays(r).close_id = id;
        rays(r).close_pt = close_pt;
        rays(r).min_dist = di;

        close_rays = [close_rays, rays(r)];

        ids = [ids, id];
    end
    if di < min_di
        min_di = di;
    end
end
if isempty(close_rays)
    minima = NaN;
    close_ray_mat = NaN;

    return
end
% fprintf('Located %i potentially valid rays, down to %1.1f km\n', length(close_rays), min_di / 1E3)


%% set up the matrix 
close_ray_mat = ones(length(elvarr), length(azarr)) * 1E9;
for r = 1:length(close_rays)
    close_ray_mat(...
        elvarr == close_rays(r).initial_elev, ...
        azarr == close_rays(r).initial_bearing) = close_rays(r).min_dist;
end


%% find minima
[az2d, el2d] = meshgrid(azarr, elvarr);

TF = imregionalmax(-close_ray_mat);
minaz = az2d(TF);
minel = el2d(TF);

for i = 1:length(minaz)
    for r = 1:length(close_rays)
        if close_rays(r).initial_elev == minel(i) && ...
            close_rays(r).initial_bearing == minaz(i)
            minima(i) = close_rays(r);
        end
    end
end

% %% plotting
% hold on
% 
% close_ray_mat(close_ray_mat == 1E9) = close_ray_mat(close_ray_mat == 1E9) * NaN;
% contourf(elvarr, azarr, close_ray_mat');
% plot(minel, minaz, 'r.', 'MarkerSize', 15)
% 
% hold off

