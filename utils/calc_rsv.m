function rsv = calc_rsv(origin_lat, origin_long, origin_ht, elv, az, freq, ...
    OX_mode, Ne, Bx_i, By_i, Bz_i)

%% calculate the ray state vector, based on the origin and bearing plus the refractive index
% rsv = calc_rsv(origin_lat, origin_long, origin_ht, elv, az, freq, ...
%           OX_mode, Ne, Bx_i, By_i, Bz_i)


%% test inputs
assert(~isnan(Ne), 'Ne should not be NaN')
%% ray geom
[rsv.pos_x, rsv.pos_y, rsv.pos_z] = wgs84_llh2xyz(origin_lat, origin_long, origin_ht * 1E3);
[rsv.dir_x, rsv.dir_y, rsv.dir_z] = relaz2xyz(1, elv, az, origin_lat, origin_long);
rsv.group_path = 0;
rsv.geom_path = 0;
rsv.phase_path = 0;
rsv.absorption = 0;
rsv.indep_var = 0;
rsv.ODE_step_size = 1000;

%% Get plasmafreq at the transmitter

plasmafreq = sqrt(80.6 * Ne ./ 1E6);

% Local plasma > transmission freq - skip
if plasmafreq > freq
    % fprintf('local plasma environment > transmit freq (%1.1e > %1.1e MHz - skipping\n', ...
    %     plasmafreq, freq)
    rsv = NaN;
    return
end

%% Determine refractive index
theta = vector_angle([rsv.dir_x, rsv.dir_y, rsv.dir_z], [Bx_i, By_i, Bz_i]);
B = sqrt(Bx_i.^2 + By_i.^2 + Bz_i.^2);
[n_O, n_X, n] = appleton_hartree(deg2rad(theta), Ne * 1E6, B, freq * 1E6);

% [pick the right one]
switch OX_mode
    case 1
        refractive_ind = n_O;
    case -1
        refractive_ind = n_X;
    case 0
        refractive_ind = n;
end

if refractive_ind > 1
    rsv = NaN;
    return
end

%% scale rsv.dir_x, rsv.dir_y, rsv.dir_z by the refractive index
rsv.dir_x = rsv.dir_x * refractive_ind;
rsv.dir_y = rsv.dir_y * refractive_ind;
rsv.dir_z = rsv.dir_z * refractive_ind;





