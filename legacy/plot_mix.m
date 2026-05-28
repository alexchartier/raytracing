%% specify inputs
Times = datenum(2019, 3, 2, 17, 40, 0):20/60/24:datenum(2019, 3, 3);
pot_fname_fmt = '/Users/chartat1/pymix/data/pot_sami_cond/mar19/ampere_mix_{YYYY}-{mm}-{DD}T{HH}-{MM}-{SS}Z.nc';
plot_fname_fmt = '/Users/chartat1/raytracing/plots/potential/ampere_mix_{YYYY}-{mm}-{DD}T{HH}-{MM}-{SS}Z.png';

lat_cutoff = 45;


%% loop over times
for t = 1%:length(Times)
    time = Times(t);

    %% Load
    pot_fname = filename(pot_fname_fmt, time);
    mix = load_mix(pot_fname);

    magpole = [unique(mix.Geographic_Latitude(:, 1)), unique(mix.Geographic_Longitude(:, 1))];

    %% Plote
    figure('Position', [0, 0, 800, 1000])
    
    hold on
    m_proj('Stereographic', 'latitude', 90, 'radius', 90 - lat_cutoff, 'rotation', 90)
    lat = [mix.Geographic_Latitude; mix.Geographic_Latitude(1, :)];
    lon = [mix.Geographic_Longitude; mix.Geographic_Longitude(1, :)];
    pot = [mix.Potential; mix.Potential(1, :)];
    [C, hC] = m_contourf(lon, lat, pot);
    clabel(C, hC)
    a = colorbar();
    colormap(bluewhitered)
    %clim([-100, 100])
    ylabel(a,'Potential (kV)','FontSize',20,'Rotation',90);
    % set(hC, 'LineStyle', 'None')

    m_grid('color', 'k', 'FontSize', 20)
    m_coast('color', 'k');
    title(filename('{YYYY/mm/dd HH:MM} UT', time))
    
    % Sun's position
    sun_offset = 0.001;
    sun_lon = - (time - floor(time)) * 360 + 180;
    sun_lon(sun_lon > 180) = sun_lon - 360;
    m_plot(sun_lon, lat_cutoff + sun_offset, '.k', 'markersize', 50)

    set(gca, 'FontSize', 20)
    hold off
    
    %% Save Plots as Files
%     Name = filename(plot_fname_fmt, time);
%     export_fig(Name);
%     
%     fprintf('Saved to %s\n', Name)
%     close


end

