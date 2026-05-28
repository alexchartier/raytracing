function tdoa_tbl = compute_tdoa_meas(t_sec, sv1_loc_x_y_z, sv2_loc_x_y_z, ground_loc_x_y_z, freq_hz, iono_map_file, varargin)
%% compute_tdoa_meas.m
% Compute time-difference-of-arrival (TDoA) using ionospheric propagation
% model, SV positional data (TLE/ECEF), and ground station position (ECEF)
%
% Required Arguments:
%   sv1_loc_z_y_z       - Space vehicle 1 coordinates (x,y,z in meters, ECEF)
%   sv2_loc_z_y_z       - Space vehicle 2 coordinates (x,y,z in meters, ECEF)
%   ground_loc_x_y_z    - Ground station coordinates (x,y,z in meters, ECEF)
%   freq_hz             - Center frequency (Hz)
%   iono_map_file       - Ionosphere map file (.dat)
% Optional Arguments (varargin):
%   missing_ray_opt     - What to do when no rays are found ('drop', 'repeat', 'interp')
%                           'drop'      - drop time point (default)
%                           'repeat'    - repeat closest value point
%                           'interp'    - interpolate between points
%   parallel            - Number of processes to use to speed up calculations ('parfor', default=1)
%   verbose             - Verbose outputs
%
% Outputs:
%   tdoa_tbl    - 'timetable' containing TDoA measurements/covariance
%                   | time (sec) | tdoa (us) | covar     |
%                   | 1          | tdoa(1)   | covar(1)  |
%                   | 2          | tdoa(2)   | covar(2)  |
%                   | ...        | ...       | ...       |
%                   | N          | tdoa(N)   | covar(N)  |
%
%                 Table may contain temporal discontinuities when using
%                 the 'drop' option for missing rays

missing_ray_opt = 'drop';
parallel        = 1;
verbose         = 0;

if length(varargin) >= 1
    missing_ray_opt = varargin{1};
end
if length(varargin) >= 2
    parallel = varargin{2};
end
if length(varargin) >= 3
    verbose = varargin{3};
end

%% Prepare
size_max = min([length(sv1_loc_x_y_z) length(sv2_loc_x_y_z)]);
dt       = abs(t_sec(2) - t_sec(1));
t_min    = t_sec(1);
t_max    = t_sec(size_max);

t = seconds((t_min:dt:t_max).');
% CubeSat1
receiver1_loc_x_y_z = sv1_loc_x_y_z; % m
% CubeSat2
receiver2_loc_x_y_z = sv2_loc_x_y_z; % m
% Calculate delta T from ionosphere (Alex code) (us)
freq_carr_mhz = freq_hz / 1e6;

%% Compute
TDOA = zeros(length(t), 1);
n    = zeros(length(t), 1);

% Vectorize ionospheric delay computation
sv_iono_delay = zeros(length(t), 2);
if parallel > 1
    parfor (i = 1:length(t), parallel)
        delay_us = calc_ionospheric_delay(iono_map_file, ground_loc_x_y_z, receiver1_loc_x_y_z(i,:), freq_carr_mhz, verbose);
        sv_iono_delay(i, 1) = delay_us(1);
    end
    parfor (i = 1:length(t), parallel)
        delay_us = calc_ionospheric_delay(iono_map_file, ground_loc_x_y_z, receiver2_loc_x_y_z(i,:), freq_carr_mhz, verbose);
        sv_iono_delay(i, 2) = delay_us(1);
    end
else
    for i = 1:length(t)
        if verbose
            fprintf("tdoa: %d/%d\n", i, length(t));
        end
        delay1_us = calc_ionospheric_delay(iono_map_file, ground_loc_x_y_z, receiver1_loc_x_y_z(i,:), freq_carr_mhz, verbose);
        delay2_us = calc_ionospheric_delay(iono_map_file, ground_loc_x_y_z, receiver2_loc_x_y_z(i,:), freq_carr_mhz, verbose);
        sv_iono_delay(i, 1) = delay1_us(1);
        sv_iono_delay(i, 2) = delay2_us(1);
    end
end

% Remove NaNs
nan_idx = sum(isnan(sv_iono_delay), 2);

sv_iono_delay = sv_iono_delay(~nan_idx,:);
t_valid       = t(~nan_idx);

% Compute TDoAs
tdoa = sv_iono_delay(:,1) - sv_iono_delay(:,2);

% Return TDoA measurements (time-aligned)
tdoa_tbl = timetable(t_valid, tdoa, VariableNames={'tdoa'});
