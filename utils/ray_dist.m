function [dl, closest_pt, id, group_path, geom_path, total_absorption] = ray_dist(ray, loc)
%% Calculate distance of ray from point
% ray - structure from raytrace_3d
% loc - lat (deg), lon (deg), height (km)
% dl - distance from point to line
% id - index of closest point on raypath
% reflect - if true, reject the pre-reflection part of the ray (should be
% false for transionospheric rays)

%% Work out if we're reflecting or not
reflect = test_reflect(ray.height(1), loc(3));

%% clean out the NAN part of the ray
finind = ~isnan(ray.height) & ray.height < 1E40;
for vn = {'lat', 'lon', 'height', 'group_range'}
    ray.(vn{1}) = ray.(vn{1})(finind);
end


%% Assess distance from ray to RX location
dl = 1E6;  % Starting ("bad") distance in case it doesn't home

[x, y, z] = wgs84_llh2xyz(ray.lat, ray.lon, ray.height * 1E3);
ray_xyz = [x; y; z];

closest_pt = nan;

id = NaN;
group_path = NaN;
geom_path = NaN;
total_absorption = NaN;

%% Optional reflection testing
if reflect
    refloc = find(diff(diff(ray.height) >= 0, 1)); % Locate post-reflection part of ray

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
    if (loc(3) - max(ray.height(refloc:end))) > 0
        %disp('Ray did not get back to receiver height')
        return
    end
    ray_xyz = ray_xyz(:, refloc:end);  % skipping pre-reflection points

else
    refloc = 1;
end


%% find closest point along ray - interpolate between closest few of them
[x, y, z] = wgs84_llh2xyz(loc(1), loc(2), loc(3) * 1E3);
loc_xyz = [x, y, z]' .* ones(3, size(ray_xyz, 2));

dists = sqrt(sum((loc_xyz - ray_xyz) .^2));

id = find(dists == min(dists));

if id == 1
    % start of ray is closest - return out
    %disp('start pt')
    dl = NaN;
    return
end

if length(ray.lat) == 2
    % weird short ray - return out
    %disp('start pt')
    dl = NaN;
    return
end

closest_pt = ray_xyz(:, id);
% if id == length(dists)
%     % end of ray is closest - return that
%     %disp('endpt')
%     dl = dists(id);
%     group_path = NaN;
%     return
% end

% select a point before closest point
id0 = id - 1;

dl = point_to_line(loc_xyz(:, 1)', ray_xyz(:, id0)', ray_xyz(:, id)');

if dl == 1E6
    dl = NaN;
    group_path = NaN;
    total_absorption = NaN;
    geom_path = NaN;
end


%% check the point really is close (point-to-line looks at an infinite line)
if sqrt(sum((closest_pt - loc_xyz(:, 1)).^2)) > 1E6
    dl = NaN;
    group_path = NaN;
    total_absorption = NaN;
    geom_path = NaN;
else

    %% find the group range to the receiver
    idx = id-2;

    if idx < 1
        idx = 1;
    end
    group_path = interp1(sqrt(sum(ray_xyz(:, idx:end).^2)), ...
        ray.group_range(idx -1 + refloc:end), sqrt(sum(loc_xyz(:, 1) .^2)),...
        'linear', 'extrap');
    
    ray.geometric_distance = ray.geometric_distance(1:length(ray.group_range));
    ray.absorption = ray.absorption(1:length(ray.group_range));

    geom_path = interp1(sqrt(sum(ray_xyz(:, idx:end).^2)), ...
        ray.geometric_distance(idx -1 + refloc:end), sqrt(sum(loc_xyz(:, 1) .^2)),...
        'linear', 'extrap');
    total_absorption = interp1(sqrt(sum(ray_xyz(:, idx:end).^2)), ...
        ray.absorption(idx -1 + refloc:end), sqrt(sum(loc_xyz(:, 1) .^2)),...
        'linear', 'extrap');


end



