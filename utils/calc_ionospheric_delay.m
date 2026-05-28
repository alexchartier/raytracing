function [ionosphere_delay_us, homed_rays] = calc_ionospheric_delay(...
    ionosphere_map_fname, ground_loc_x_y_z, receiver_loc_x_y_z, freq_carr_mhz)
%% calc_ionospheric_delay.m
% Calculates ionospheric delays from model ionosphere and tx/rx locations
% Note this delay includes both ionospheric and non-ionospheric (geometric)
% contributions. 
%
% ionosphere_map_fname - Ionospheric model file in SAMI3 format
% ground_loc_x_y_z - transmitter location in ECEF (m)
% receiver_loc_x_y_z - receiver location in ECEF (m)
% freq_carr_mhz - Transmitted signal frequency
% 
% ionosphere_delay_us - array of delays in microseconds (n_OX_mode X n_paths) 
% OX_mode defined as [-1, 1] (X and O)
% 
% Example:
%
% time = datenum(2017, 1, 10, 18, 0, 0);
% in_fn_fmt = '/Users/chartat1/data/sami3/2017_tid/sami_mat/{YYYY-mm-dd_HHMM}.mat';
% ionosphere_map_fname = filename(in_fn_fmt, time);
% RE = 6371E3;
% ground_loc_x_y_z = sphcart([RE, deg2rad(38), deg2rad(283)]);
% receiver_loc_x_y_z = sphcart([RE + 400E3, deg2rad(35), deg2rad(283)]);
% freq_carr_mhz = 10;
% ionosphere_delay_us = calc_ionospheric_delay(ionosphere_map_fname, ...
%     ground_loc_x_y_z, receiver_loc_x_y_z, freq_carr_mhz)

%% Load/inputs
maxdist = 1E5;  % meters from homing
tol = [1e-7 0.01 25];       % ODE solver tolerance and min max stepsizes
homing_tol_m = 750;
OX_modes = [-1, 1];
c = 299792458;

txloc = convert_coords(ground_loc_x_y_z);
rxloc = convert_coords(receiver_loc_x_y_z);

fprintf("Transmitter: %1.1f °N %1.1f °E %1.1f km\n", txloc(1), txloc(2), txloc(3))
fprintf("Receiver:    %1.1f °N %1.1f °E %1.1f km\n", rxloc(1), rxloc(2), rxloc(3))

%%
model = loadstruct(ionosphere_map_fname);
% [iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
%     Bx, By, Bz, geomag_grid_parms] = gen_grid_parms(model);

%% Raytrace

for mi = 1:length(OX_modes)
    rays = gen_ionogram(freq_carr_mhz, OX_modes(mi), txloc, rxloc, ...
        model.iono_en_grid, model.iono_en_grid, model.collision_freq, ...
        model.iono_grid_parms, model.geomag_grid_parms, ...
        model.Bx, model.By, model.Bz, ...
        maxdist, tol, homing_tol_m);
    if length(rays) > 0
    homed_rays{mi} = rays{1};
    end
end

%% Store out the delays
maxlen = 0;

if exist('homed_rays', 'Var') == 1
    for i = 1:length(homed_rays)
        if length(homed_rays{i}) > maxlen
            maxlen = length(homed_rays{i});
        end
    end
    ionosphere_delay_us = zeros(length(OX_modes), maxlen) * NaN;
    for i = 1:length(homed_rays)
        for ri = 1:length(homed_rays{i})
            ionosphere_delay_us(i, ri) = homed_rays{i}(ri).group_range_to_rx * 1E3 / c * 1E6;
        end
    end

else
    ionosphere_delay_us = NaN;
    homed_rays = NaN;
end


















