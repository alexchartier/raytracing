%% plot_twosat_osse.m
% script to plot the output of twosat_synth/sat_synth_nebula

%% Paths

%datadir = '~/data/nebula/iri_osse_el80-0_az0-360_2sats_hilat/homings/sat1_sat2/';
datadir = '~/data/nebula/iri_osse_el80-0_az0-360_2sats_lowlat/homings/sat1_sat2/';
fontsize = 24;

%% plots
flist = dir(datadir);
close all

for fi = 3% :length(flist)
    fn = [flist(fi).folder, '/', flist(fi).name];
    rays = loadstruct(fn);
    OX_mode = rays{1}(1).OX_mode;
    switch OX_mode
        case -1
            OX_mode_name = 'X';
        case 1
            OX_mode_name = 'O';
        case 0
            OX_mode_name = 'No_B';
    end

    lat = mean(rays{end}(1).lat);
    lon = mean(rays{end}(1).lon);
    fmin = rays{1}(1).frequency;

    title_str = sprintf('%s mode, %1.1f N, %1.1f E, Fmin: %1.1f MHz, Spacing: 600 km', ...
        OX_mode_name, lat, lon, fmin);
    
    plot_ionogram(rays, title_str)
    set(gca, 'FontSize', fontsize, 'FontName', 'Futura')
    figure
    plot_rays_on_2d_ionosphere(rays, '');

end

%% minimum freq and maximum range
fmin = 30;
rmax = 0;
for fi = 3 :length(flist)
    fn = [flist(fi).folder, '/', flist(fi).name];
    rays = loadstruct(fn);
    if rays{1}(1).frequency < fmin
        fmin = rays{1}(1).frequency;
    end
    for ri = 1:length(rays)
        for r = 1:length(rays{ri})
            if rays{ri}(r).group_range_to_rx > rmax 
                rmax = rays{ri}(r).group_range_to_rx;
            end
        end
    end

end

fprintf('Fmin = %1.1f MHz\nRmax = %1.1f km\n', fmin, rmax)




