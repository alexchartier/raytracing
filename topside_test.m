%% raytrace_3d test from topside
UT = [2000 9 21 0 0];           % UT - year, month, day, hour, minute
speed_of_light = 2.99792458e8;
R12 = 100;
elvarr = -90:3:-30;               % initial elevation of rays
freq = 8;   % frequency (MHz)
azarr = -180:3:180; % initial bearing of rays
origin_lat = -10.0;             % latitude of the start point of rays
origin_long = 130.0;            % longitude of the start point of rays
origin_ht = 800.0;              % altitude of the start point of rays
doppler_flag = 1;               % interested in Doppler shift


[elv2d,az2d] = meshgrid(elvarr,azarr);
elv1d = double(elv2d(:))';
az1d = double(az2d(:))';
freqs = ones(size(elv1d))*freq;

fprintf( ['\n' ...
   'Topside example of 3D magneto-ionic numerical raytracing for a WGS84 ellipsoidal' ...
   ' Earth\n\n'])


%% generate ionospheric, geomagnetic and irregularity grids
ht_start = 100;          % start height for ionospheric grid (km)
ht_inc = 2;             % height increment (km)
num_ht = 401;           
lat_start = -20.0;
lat_inc = 0.3;
num_lat = 101.0;
lon_start= 128.0;
lon_inc = 1.0;
num_lon = 5.0;
iono_grid_parms = [lat_start, lat_inc, num_lat, lon_start, lon_inc, num_lon, ...
      ht_start, ht_inc, num_ht, ];

B_ht_start = ht_start;          % start height for geomagnetic grid (km)
B_ht_inc = 10;                  % height increment (km)
B_num_ht = ceil(num_ht .* ht_inc ./ B_ht_inc);
B_lat_start = lat_start;
B_lat_inc = 1.0;
B_num_lat = ceil(num_lat .* lat_inc ./ B_lat_inc);
B_lon_start = lon_start;
B_lon_inc = 1.0;
B_num_lon = ceil(num_lon .* lon_inc ./ B_lon_inc); 
geomag_grid_parms = [B_lat_start, B_lat_inc, B_num_lat, B_lon_start, ...
      B_lon_inc, B_num_lon, B_ht_start, B_ht_inc, B_num_ht];

tic
fprintf('Generating ionospheric and geomag grids... ')
[iono_pf_grid, iono_pf_grid_5, collision_freq, Bx, By, Bz] = ...
    gen_iono_grid_3d(UT, R12, iono_grid_parms, ...
                     geomag_grid_parms, doppler_flag);
toc
fprintf('\n')

% convert plasma frequency grid to  electron density in electrons/cm^3
iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;

%
%% call raytrace
%
nhops = 1;                  % number of hops
tol = [1e-7 0.01 25];       % ODE solver tolerance and min max stepsizes
num_elevs = length(elv1d);

% Create starting raypath state vector  
clear ray_state_vec_in
    for i = 1:length(elv1d)
        [rsv.pos_x, rsv.pos_y, rsv.pos_z] = wgs84_llh2xyz(origin_lat, origin_long, origin_ht * 1E3);
        [rsv.dir_x, rsv.dir_y, rsv.dir_z] = relaz2xyz(1, elv1d(i), az1d(i), origin_lat, origin_long);
        rsv.group_path = 0;
        rsv.geom_path = 0;
        rsv.phase_path = 0;
        rsv.absorption = 0;
        rsv.indep_var = 0;
        rsv.ODE_step_size = 1000;
        ray_state_vec_in(i) = rsv;
    end
    
    

% Generate the O mode rays
OX_mode = 2;
fprintf('Generating %d O-mode rays ...', num_elevs);
tic
[ray_data_O, ray_O, ray_state_vec_O] = ...
  raytrace_3d(origin_lat, origin_long, origin_ht, elv1d, az1d, freqs, ...
              OX_mode, nhops, tol, iono_en_grid, iono_en_grid_5, ...
	          collision_freq, iono_grid_parms, Bx, By, Bz, ...
	          geomag_grid_parms, ray_state_vec_in);
	      
NRT_total_time = toc;
fprintf('\n   NRT-only execution time = %f, Total mex execution time = %f\n\n', ...
        [ray_data_O.NRT_elapsed_time], NRT_total_time)

    
%% Plot

% plot the rays
figure(1)
hold on
pos = get(gcf, 'position');
pos(3) = pos(3)*1.5;
pos(4) = pos(4)*1.5;
set(gcf, 'position', pos)
plot3(ray_O(1).lat, mod(ray_O(1).lon, 360), ray_O(1).height, '.b', ...
      'markersize', 5)
set(gca, 'Zlim', [150 850])

for ii = 1:num_elevs
  plot3(ray_O(ii).lat, mod(ray_O(ii).lon, 360), ray_O(ii).height, '.b', ...
        'markersize', 5)
end
hold off

