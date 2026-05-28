%% specify inputs
Times = datenum(2019, 3, 2, 17, 40, 0):20/60/24:datenum(2019, 3, 3);
mod_fname_fmt = '/Users/chartat1/data/sami3/2019/sami3_regulargrid_ancillary_{YYYYmmmDD}.nc';
plot_fname_fmt = '/Users/chartat1/raytracing/plots/nmf2/nmf2_{YYYY}-{mm}-{DD}T{HH}-{MM}-{SS}Z.png';

lat_cutoff = 45;
vars = {'lat', 'lon', 'nmf2', 'time'};



%% loop over times
for t = 1%:length(Times)
    time = Times(t);

    %% Load
    if time == floor(time) || time == Times(1)
        mod_fname = filename(mod_fname_fmt, time);
        sami = load_sami(mod_fname, vars);
    end
    sami.time = datenum(sami.time);
        
    sami_time = closest(sami.time, time);
    tind = sami.time == sami_time;
    sami_t = sami;
    sami_t.nmf2 = sami.nmf2(:, :, tind);

    %% Plot
    figure('Position', [0, 0, 800, 1000])
    hold on
    m_proj('Stereographic', 'latitude', 90, 'radius', 90 - lat_cutoff, 'rotation', 90)
    lat = [sami_t.lat];
    lon = [sami_t.lon; sami_t.lon(1)];
    nmf2 = [sami_t.nmf2; sami_t.nmf2(1, :)] / 1E5;
    [C, hC] = m_contourf(lon, lat, nmf2');
    clabel(C, hC)
    a = colorbar();
    clim([0, 5])
    ylabel(a,'NmF2 (1E5 el. m-3)','FontSize',20,'Rotation',90);
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
    Name = filename(plot_fname_fmt, time);
    export_fig(Name);
    
    fprintf('Saved to %s\n', Name)
    close


end

