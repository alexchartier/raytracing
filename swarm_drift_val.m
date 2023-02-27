%% swarm_drift_val.m
% Validate the MIX potential against Swarm ExB drift

%% specify inputs
Times = datenum(2019, 3, 2, 10, 0, 0):2/60/24:datenum(2019, 3, 3);
    
mlat_cutoff = 60;
pot_fname_fmt = '/Users/chartat1/pymix/data/pot_sami_cond/mar19/ampere_mix_{YYYY}-{mm}-{DD}T{HH}-{MM}-{SS}Z.nc';
plot_fname_fmt = '/Users/chartat1/pymix/plots/potential/ampere_mix_{YYYY}-{mm}-{DD}T{HH}-{MM}-{SS}Z.png';

swarm_fn = '~/raytracing/data/swarm_proc/2019Mar02_swarm.nc';

sats = ['A', 'B'];
dec_rate =  10;

ref_rad = 6800E3;
lat_cutoff = 40;

%% Loop over time
time = Times(20);

%% Load MIX potential 
close
pot_fname = filename(pot_fname_fmt, time);
mix = [];

keys = {'Potential', 'Geographic Latitude', 'Geographic Longitude', 'MLAT (AACGM)', 'MLON (AACGM)'};

for k = keys
    outname = regexprep(k{1}, ' ', '_');
    outname = regexprep(outname, '(', '');
    outname = regexprep(outname, ')', '');
    mix.(outname) = double(ncread(pot_fname, k{1}));
    mix.(outname) = mix.(outname)(2:end, :); 
end
[mix.efld_N_Mag, mix.efld_E_Mag] = spherical_gradient(...
    ref_rad, mix.MLAT_AACGM, mix.MLON_AACGM, mix.Potential);

mix.efld_N_Mag = mix.efld_N_Mag .* 2 .* pi .* ref_rad; 
mix.efld_E_Mag = mix.efld_E_Mag .* 2 .* pi .* ref_rad .* cosd(mix.MLAT_AACGM);

magpole = [unique(mix.Geographic_Latitude(:, 1)), unique(mix.Geographic_Longitude(:, 1))];


%% Load Swarm ExB drift

% ncdisp(swarm_fn)

keys = {'Latitude', 'Longitude', 'Radius', 'Vn', 'Ve', 'POSIXtime'};
for k = keys
    swarm.(k{1}) = ncread(swarm_fn, k{1});
    swarm.(k{1}) = swarm.(k{1})(1:dec_rate:end);
end
swarm.time = datetime(swarm.POSIXtime, 'ConvertFrom', 'posixtime' );    

%%
tidx = within(datenum(swarm.time), time, 5/60/24);
hold on
m_proj('Stereographic', 'latitude', 90, 'radius', 90 - lat_cutoff, 'rotation', 90)
[~, hC] = m_contourf(lon, lat, pot');
colorbar()
% set(hC, 'LineStyle', 'None')

m_quiver(swarm.Longitude(tidx),swarm.Latitude(tidx),swarm.Ve(tidx), swarm.Vn(tidx),10, 'k');
m_grid('color', 'k', 'FontSize', 20)
title(filename('{YYYY mm dd HH:MM}', time))
set(gca, 'FontSize', 20)
hold off
%% Save Plots as Files
% Name = filename(OPath,Time);
% export_fig(Name);


% figure
% m_proj('Stereographic', 'latitude', 90, 'radius', 90 - lat_cutoff, 'rotation', 90)
% m_contourf(lon_out, lat_out, pot_out')
% m_grid('color', 'k', 'FontSize', 14)


%% Calculate ExB drifts using spherical gradient


%% Compare against Swarm cross-track ExB drift
    












