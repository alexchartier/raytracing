%% Prepare
close all;
clear;

%% Inputs 
% SV TLE files (position vs time, ECEF format)
sc_crd_fn1 = '~/Downloads/CubeSat1_Fixed_Position_Velocity.txt';
sc_crd_fn2 = '~/Downloads/CubeSat2_Fixed_Position_Velocity.txt';

% Ionospheric data
ionosphere_map_fname = 'files/2017-01-10_1800.mat';

% Initial ground location of known emitter (ECEF format)
% NOTE: Comment out the following variable to skip the initial TDoA step
ground_loc_x_y_z = [1089.221557, -4913.160938, 3905.460202] * 1e3; % x,y,z meters

% Geolocated ground coordinates, from geo engine
% NOTE: Comment out the following variable to skip the estimated TDoA step
% Iteration 1 (working)
% geo_ground_loc_x_y_z = [1213547.7322280437,-4943043.3654513685, 3921976.831505978];
% Iteration 2 (no valid rays)
geo_ground_loc_x_y_z = [1068512.5857582283, -4905159.662085281, 3906876.5075383973];

% Carrier frequency
freq_carr_mhz = 10;

% Valid time slice from file (indices)
% NOTE: Valid window from ~9-11 minutes (empirically determined)
t_sec = 9*60 : 11*60;

% Computation parallelism (1 for single thread)
parallel = 4;

% Use cached results (a.k.a. use files if they exist)
use_cache = 0;

% Plot results
use_plots = 1;



%% Read files
% SV TLE files
CubeSat1_file = readmatrix(sc_crd_fn1);
CubeSat2_file = readmatrix(sc_crd_fn2);

% CubeSat1
CubeSat1.x_pos = CubeSat1_file(:,5) * 1000; % m
CubeSat1.y_pos = CubeSat1_file(:,6) * 1000; % m
CubeSat1.z_pos = CubeSat1_file(:,7) * 1000; % m
receiver1_loc_x_y_z = [CubeSat1.x_pos, CubeSat1.y_pos, CubeSat1.z_pos]; % m

% CubeSat2
CubeSat2.x_pos = CubeSat2_file(:,5) * 1000; % m
CubeSat2.y_pos = CubeSat2_file(:,6) * 1000; % m
CubeSat2.z_pos = CubeSat2_file(:,7) * 1000; % m
receiver2_loc_x_y_z = [CubeSat2.x_pos, CubeSat2.y_pos, CubeSat2.z_pos]; % m

% SV positions w/r/t valid time slice
receiver1_loc_x_y_z = receiver1_loc_x_y_z(t_sec,:);
receiver2_loc_x_y_z = receiver2_loc_x_y_z(t_sec,:);

% Output file name prefix
tdoa_file_prefix = sprintf("tdoa_meas_profile_2_%.3fmhz", freq_carr_mhz);


%% Compute initial TDoAs (ionospheric model)
% use initial coordinates

if exist('ground_loc_x_y_z', 'var')
    tdoa_fname = sprintf("%s_raw.csv", tdoa_file_prefix);
    
    if use_cache && exist(tdoa_fname, 'file')
        % Read TDoAs from file
        fprintf("reading ionospheric delay (cache \'%s\')...\n", tdoa_fname);
        tdoa_tbl = readtimetable(tdoa_fname);
    else
        % Compute TDoAs w/ ionospheric delay
        fprintf("computing ionospheric delay...\n");
        tic
        tdoa_tbl = compute_tdoa_meas(t_sec, receiver1_loc_x_y_z, receiver2_loc_x_y_z, ground_loc_x_y_z, freq_carr_mhz*1e6, ionosphere_map_fname, 'drop', parallel, 0);
        elapsed = toc;
        fprintf("%d TDoAs in %.3f secs (%.3f TDoAs/sec)\n", length(t_sec), elapsed, length(t_sec)/elapsed);
    
        % Store to file
        tdoa_tbl.Properties.RowTimes.Format = 'hh:mm:ss.SSS';
        writetimetable(tdoa_tbl, tdoa_fname);
    end
    
    if isempty(tdoa_tbl)
        fprintf("[WARNING] initial TDoA vector empty (no rays found), please check input parameters\n");
    end

        % Compute geometric TDoAs (free-space)
    fprintf("computing free-space path delay...\n");
    
    % compensated geometric TDoAs (w/ geo engine position)
    tdoa_geom_tbl = compute_tdoa_geom(t_sec, receiver1_loc_x_y_z, receiver2_loc_x_y_z, ground_loc_x_y_z, 1);
    
    if use_plots
        figure();
        plot(tdoa_tbl.Properties.RowTimes, tdoa_tbl.tdoa, 'b-o');
        hold on;
        plot(tdoa_geom_tbl.Properties.RowTimes, tdoa_geom_tbl.tdoa, 'k-.');
        grid on;
        title_str = sprintf("Emitter @ %.3e,%.3e,%.3e", ...
            ground_loc_x_y_z(1), ground_loc_x_y_z(2), ground_loc_x_y_z(3));
        title("Initial TDoAs");
        subtitle(title_str);
        xlabel("Time [sec]");
        ylabel("TDoA (SV_1 - SV_2) [us]");
    end
end


%% Call geolocation engine w/ initial TDoAs
% geo_ground_loc_x_y_z = call_geo_tool(tdoa_fname);



%% Compute updated TDoAs (ionospheric model)
% use geolocated ground coordinates (from geo engine)

if exist('geo_ground_loc_x_y_z', 'var')
    tdoa_fname = sprintf("%s_comp.csv", tdoa_file_prefix);
    
    if use_cache && exist(tdoa_fname, 'file')
        fprintf("reading ionospheric delay w/ estimated coords (cache \'%s\')...\n", tdoa_fname);
        tdoa_comp_tbl = readtimetable(tdoa_fname);
    else
        fprintf("computing ionospheric delay w/ estimated coords...\n");
        tic
        tdoa_comp_tbl = compute_tdoa_meas(t_sec, receiver1_loc_x_y_z, receiver2_loc_x_y_z, geo_ground_loc_x_y_z, freq_carr_mhz*1e6, ionosphere_map_fname, 'drop', parallel, 0);
        elapsed = toc;
        fprintf("%d TDoAs in %.3f secs (%.3f TDoAs/sec)\n", length(t_sec), elapsed, length(t_sec)/elapsed);
    
        % Store to file
        tdoa_comp_tbl.Properties.RowTimes.Format = 'hh:mm:ss.SSS';
        writetimetable(tdoa_comp_tbl, tdoa_fname);
    end
    
    if isempty(tdoa_comp_tbl)
        fprintf("[WARNING] compensated TDoA vector empty (no rays found), please check input parameters\n");
    end
    
    
    % Compute geometric TDoAs (free-space)
    fprintf("computing free-space path delay w/ estimated coords...\n");
    
    % compensated geometric TDoAs (w/ geo engine position)
    tdoa_geom_comp_tbl = compute_tdoa_geom(t_sec, receiver1_loc_x_y_z, receiver2_loc_x_y_z, geo_ground_loc_x_y_z, 1);
    
    % Store to file
    tdoa_fname = sprintf("%s_comp_geom.csv", tdoa_file_prefix);
    tdoa_geom_comp_tbl.Properties.RowTimes.Format = 'hh:mm:ss.SSS';
    writetimetable(tdoa_geom_comp_tbl, tdoa_fname);
    
    
    % Compute normalized TDoA measurements
    
    fprintf("computing TDoAs w/ estimated coords...\n");
    % Remove missing rays from results ('NaN' in data)
    tdoa_cmp_tbl   = rmmissing(synchronize(tdoa_comp_tbl, tdoa_geom_comp_tbl));
    delta_tdoa     = tdoa_cmp_tbl.tdoa_tdoa_comp_tbl - tdoa_cmp_tbl.tdoa_tdoa_geom_comp_tbl;
    tdoa_final_tbl = timetable(tdoa_cmp_tbl.Properties.RowTimes, ...
                               tdoa_cmp_tbl.tdoa_tdoa_comp_tbl, ...
                               tdoa_cmp_tbl.tdoa_tdoa_geom_comp_tbl, ...
                               delta_tdoa, ...
                               VariableNames={'tdoa', 'tdoa_comp', 'delta_tdoa'});

    % Store to file
    tdoa_fname = sprintf("%s_delta.csv", tdoa_file_prefix);
    tdoa_final_tbl.Properties.RowTimes.Format = 'hh:mm:ss.SSS';
    writetimetable(tdoa_final_tbl, tdoa_fname);
    
    if use_plots
        max_tdoa = max([max(abs(tdoa_final_tbl.tdoa)) max(abs(tdoa_final_tbl.tdoa_comp))]);

        figure();
        subplot(2,1,1);
        plot(tdoa_final_tbl.Properties.RowTimes, tdoa_final_tbl.tdoa,      'g-o');
        hold on;
        plot(tdoa_final_tbl.Properties.RowTimes, tdoa_final_tbl.tdoa_comp, 'k-.');
        grid on;
        title_str = sprintf("Emitter @ %.3e,%.3e,%.3e", ...
            geo_ground_loc_x_y_z(1), geo_ground_loc_x_y_z(2), geo_ground_loc_x_y_z(3));
        title("Updated TDoAs");
        subtitle(title_str);
        ylabel("TDoA_{est} (SV_1 - SV_2) [us]");
        ylim([-max_tdoa, max_tdoa]);
        legend(['TDoA_{iono}'; 'TDoA_{geom}']);
    
        subplot(2,1,2);
        plot(tdoa_final_tbl.Properties.RowTimes, tdoa_final_tbl.delta_tdoa, 'r-o');
        grid on;
        subtitle("TDoA_{iono} - TDoA_{geom}");
        xlabel("Time [sec]");
        ylabel("\Delta TDoA [us]");
        ylim([-max_tdoa, max_tdoa]);
    end

end


%% Compute compensated TDoAs (ionospheric model)
% use geolocated ground coordinates (from geo engine)

tdoa_init_fname = sprintf("%s_raw.csv",   tdoa_file_prefix);
tdoa_delt_fname = sprintf("%s_delta.csv", tdoa_file_prefix);

fprintf("computing compensated TDoAs...\n");

tdoa_init_tbl = readtimetable(tdoa_init_fname);
tdoa_delt_tbl = readtimetable(tdoa_delt_fname);
tdoa_comb_tbl = rmmissing(synchronize(tdoa_init_tbl, tdoa_delt_tbl));

tdoa_meas = tdoa_comb_tbl.tdoa_tdoa_init_tbl - tdoa_comb_tbl.delta_tdoa;
tdoa_meas_tbl = timetable(tdoa_comb_tbl.Properties.RowTimes, ...
                    tdoa_comb_tbl.tdoa_tdoa_init_tbl, ...
                    tdoa_comb_tbl.tdoa_comp, ...
                    tdoa_comb_tbl.delta_tdoa, ...
                    tdoa_meas, ...
                    VariableNames={'tdoa_init', 'tdoa_geo', 'delta_tdoa', 'tdoa_comp'});

% Store to file
tdoa_fname = sprintf("%s_final.csv", tdoa_file_prefix);
tdoa_meas_tbl.Properties.RowTimes.Format = 'hh:mm:ss.SSS';
writetimetable(tdoa_meas_tbl, tdoa_fname);
    
if use_plots
    max_tdoa = max([max(abs(tdoa_meas_tbl.tdoa_init)) max(abs(tdoa_meas_tbl.tdoa_geo))]);

    figure();
    subplot(2,1,1);
    plot(tdoa_meas_tbl.Properties.RowTimes, tdoa_meas_tbl.tdoa_init, 'b-.');
    hold on;
    plot(tdoa_meas_tbl.Properties.RowTimes, tdoa_meas_tbl.tdoa_geo,  'g-.');
    plot(tdoa_meas_tbl.Properties.RowTimes, tdoa_meas_tbl.tdoa_comp, 'k-o');
    grid on;
    title_str = sprintf("Emitter @ %.3e,%.3e,%.3e", ...
        ground_loc_x_y_z(1), ground_loc_x_y_z(2), ground_loc_x_y_z(3));
    title("TDoAs");
    subtitle(title_str);
    ylabel("TDoA (SV_1 - SV_2) [us]");
    ylim([-max_tdoa, max_tdoa]);
%     legend(['TDoA_{init}'; 'TDoA_{geo}'; 'TDOA_{comp}']);

    subplot(2,1,2);
    plot(tdoa_meas_tbl.Properties.RowTimes, tdoa_meas_tbl.delta_tdoa, 'r-o');
    grid on;
    subtitle("TDoA_{iono} - TDoA_{geom}");
    xlabel("Time [sec]");
    ylabel("\Delta TDoA [us]");
    ylim([-max_tdoa, max_tdoa]);
end