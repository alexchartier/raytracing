function [dl, closest_pt, id, group_path] = ray_dist(ray, loc, reflect)
%% Calculate distance of ray from point
% ray - structure from raytrace_3d
% loc - lat (deg), lon (deg), height (km)
% dl - distance from point to line
% id - index of closest point on raypath
% reflect - if true, reject the pre-reflection part of the ray (should be
% false for transionospheric rays)


%% Assess distance from ray to RX location
dl = 1E6;  % Starting ("bad") distance in case it doesn't home

[x, y, z] = wgs84_llh2xyz(ray.lat, ray.lon, ray.height * 1E3);
ray_xyz = [x; y; z];
% ray.x = ray_xyz(1, :);
% ray.y = ray_xyz(2, :);
% ray.z = ray_xyz(3, :);

closest_pt = nan;

id = NaN;

%% Optional reflection testing
if reflect
    refloc = find(diff(diff(ray.height) >= 0, 1));% Locate post-reflection part of ray

    %% Check if ray reflected
    if isempty(refloc)
        %disp('Ray did not reflect')
        return
    end

    if length(refloc) > 1
        refloc = refloc(1);
        %disp('multi-reflection, weird...')
    end

    %% Check if ray made it back to receiver height
    if loc(3) - max(ray.height(refloc:end)) > 0

        % disp('Ray did not get back to receiver height')
        return
    end
    ray_xyz = ray_xyz(:, refloc:end);  % skipping pre-reflection points

else
    refloc = 1;

end


%% find closest point along ray - interpolate between closest few of them
[x, y, z] = wgs84_llh2xyz(loc(1), loc(2), loc(3) * 1E3);
loc_xyz = [x, y, z]' .* ones(3, length(ray_xyz));

dists = sqrt(sum((loc_xyz - ray_xyz) .^2));

id = find(dists == min(dists));
if id == 1 
    % start of ray is closest - return out
    dl = NaN;
    return
end

% select a point before closest point
id0 = id - 1;

dl = point_to_line(loc_xyz(:, 1)', ray_xyz(:, id0)', ray_xyz(:, id)');

if dl == 1E6
    dl = NaN;
end

closest_pt = ray_xyz(:, id);
id = id + refloc  - 1;

%% find the group range to the receiver
group_path = interp1(ray.height(id-1:end), ray.group_range(id -1:end), loc(3));






























