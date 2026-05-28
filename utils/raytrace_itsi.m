function rays_out = ...
    raytrace_itsi(freq, OX_mode, txloc, rxloc, iono_en_grid, iono_en_grid_5, ...
    collision_freq, iono_grid_parms, Bx, By, Bz, geomag_grid_parms, ...
    elvarr, azarr, maxdist, tol, blind_range, homing_tol_m, varargin)

%% Raytrace_itsi
% Trace rays in the topside (may also work outside ionosphere)
%
% homed_rays = ...
%     raytrace_itsi(freq, OX_mode, txloc, rxloc, iono_en_grid, iono_en_grid_5, ...
%     collision_freq, iono_grid_parms, Bx, By, Bz, geomag_grid_parms, ...
%     elvarr, azarr, maxdist, tol, blind_range, homing_tol_m)

%% 
if length(varargin) == 1
    verbose = varargin{1};
else
    verbose = true;
end

if max(iono_en_grid) > 1E8
    disp('Electron density looks too high  - should be in el/cm3')
end


assert(txloc(3) >= 0 & rxloc(3) >= 0, 'Transmitter and receiver must be above ground')

%% PHaRLAP doesn't like to start exactly on a gridline, fix... 
[lats, lons, alts] = parms_to_lla(iono_grid_parms);  % lla of Ne

if sum(lats == txloc(1)) > 0
    txloc(1) = txloc(1) + 1E-6;
end
if sum(lons == txloc(2)) > 0
    txloc(2) = txloc(2) + 1E-6;
end

%% calculate Ne and B at the transmitter
if txloc(3) < min(alts)
    Ne = 0;
    Bx_i = 0;
    By_i = 0;
    Bz_i = 0;
else
    txloc_temp = txloc;
    if min(lons)>0 && txloc_temp(2) < 0
        txloc_temp(2) = txloc_temp(2)+360;
    end
    
    Ne = interp_3d(lats, lons, alts, iono_en_grid, txloc_temp);
    assert(~isnan(Ne), 'Interpolation failed')
    [B_lat, B_lon, B_ht] = parms_to_lla(geomag_grid_parms);  % lla of B
    if min(B_lon)>0 && txloc_temp(2) < 0
        disp('look further - switch all to txloc_temp?')
    end
    Bx_i = interp_3d(B_lat, B_lon, B_ht, Bx, txloc);
    By_i = interp_3d(B_lat, B_lon, B_ht, By, txloc);
    Bz_i = interp_3d(B_lat, B_lon, B_ht, Bz, txloc);
end


%% call raytrace to "global search" the problem
% maxdist is the maximum distance in metres
origin_lat = txloc(1);
origin_long = txloc(2);
origin_ht = txloc(3);

nhops = 1; % number of hops

plasma_freq = elec2freq(Ne)/1E3;
if plasma_freq < freq
    rays = gs_raytrace(elvarr, azarr, freq, nhops, OX_mode, ...
        origin_lat, origin_long, origin_ht, iono_en_grid, iono_en_grid_5, ...
        collision_freq, iono_grid_parms, Bx, By, Bz, geomag_grid_parms, tol, ...
        Ne, Bx_i, By_i, Bz_i);

else
    fprintf('No rays produced: local plasma freq (%1.1f) > freq (%1.1f)\n', ...
        plasma_freq, freq)
    rays_out = NaN;
    return
end

if ~isstruct(rays)
    fprintf('No rays produced for unknown reasons at %1.1f MHz\n',freq)
    rays_out = NaN;
    return
end


%% Identify the closest rays
[min_rays, ~] = cluster_rays(elvarr, azarr, rays, rxloc, maxdist);
if ~isstruct(min_rays)
    disp('No valid rays')
    rays_out = NaN;
    return
end


%% Kick out those rays with range < blind_range
good_rays = [];
for r = 1:length(min_rays)
    return_range = min_rays(r).group_range(min_rays(r).close_id);
    if return_range > (blind_range / 1E3)
        good_rays = [good_rays, min_rays(r)];
    end
end

%  fprintf('Located %i rays > blind range of %1.1f km\n', length(good_rays), blind_range/ 1E3)
if isempty(good_rays)
    disp('all rays within blind range')
    rays_out = NaN;
    return
end


%% Home from the best ray(s)
ct = 0;

elv_tol_deg = unique(diff(elvarr));
az_tol_deg = unique(diff(azarr));

for r = 1:length(good_rays)
    start_ray = rmfield(good_rays(r), {'close_id', 'close_pt'});

    homed_ray = home_ray(start_ray, rxloc, OX_mode, nhops, tol, homing_tol_m, ...
        Ne, Bx_i, By_i, Bz_i);
    
    if isstruct(homed_ray)
        ct = ct + 1;
        homed_ray.OX_mode = OX_mode;
        homed_ray.txloc = txloc;
        homed_rays(ct) = homed_ray;

    end
end

if ~exist('homed_rays', 'var')
    fprintf('Failed to home\n')
    rays_out = NaN;
    return

end

%% Go through the rays and eliminate near-duplicates
el = zeros(size(homed_rays));
az = zeros(size(homed_rays));
group_rg = zeros(size(homed_rays));

for r = 1:length(homed_rays)
    el(r) = homed_rays(r).initial_elev;
    az(r) = homed_rays(r).initial_bearing;
    group_rg(r) = homed_rays(r).group_range_to_rx;
end

first_ind_el = [];
first_ind_az = [];
first_ind_rg = [];
for i = 1:length(el)
    first_ind_el = [first_ind_el, find(round(el * 1E2) == round(el(i) * 1E2), 1)];
    first_ind_az = [first_ind_az, find(round(az * 1E2) == round(az(i) * 1E2), 1)];
    first_ind_rg = [first_ind_rg, find(round(group_rg) == round(group_rg(i)), 1)];
end



unique_id = unique(first_ind_el(ismember(first_ind_el(...
    ismember(first_ind_el, first_ind_rg)), first_ind_rg)));

rays_out = homed_rays(unique_id);

% for r = 1:length(rays_out)
% 
%     mean_gri = mean(rays_out(r).group_refractive_index);
%     max_gri =  max(rays_out(r).group_refractive_index);
%     % fprintf('mean_gri: %1.1f, max_gri: %1.1f\n', mean_gri, max_gri)
%     if mean_gri > 10
%         % disp("Something is wrong!! ")
%         rays_out(r).high_group_ref_ind = true;
%     end
% end
% 


if verbose
    fprintf('Found %i valid rays\n', length(rays_out))
end






























