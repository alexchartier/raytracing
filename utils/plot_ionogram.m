function plot_ionogram(homed_rays)
% plot_ionogram(homed_rays)

arguments 
    homed_rays cell
end


%%
colormap cool

xlimit = [0, homed_rays{end}(1).frequency];
ylimit = [0, 4000];
% climit = [0 20];
climit = [-90 0];
% climit = [0 360];
ct = 0;

hold on
freqs = [];
for ri = 1:length(homed_rays)
    rays = homed_rays{ri};
    for r = 1:length(rays)
        if ~isstruct(rays(r))
            continue
        end

        freq = rays(r).frequency;
        freqs = [freqs; freq];
        rg = rays(r).group_range_to_rx;
        if isnan(rg)
            disp(ri)
            disp(r)
            continue
        end
        % loss = rays(r).absorption(end) + fspl(rg * 1E3, 3E8 / (freq * 1E6));
        % scatter(freq, rg, 150, -loss, 'filled', 'o')

        elevation = rays(r).initial_elev; 
        azimuth = rays(r).initial_bearing; 
        absorption = rays(r).total_absorption; 

        % scatter(freq, rg, 150, absorption, 'filled', 'o')
        scatter(freq, rg, 50, elevation, 'filled', 'o')
        % scatter(freq, rg, 150, azimuth, 'filled', 'o')

        % plot(freq, rg, 'k.', 'MarkerSize', 30)
    end
end

ylim(ylimit)
xlim(xlimit) % [min(freqs(:)) - 1, max(freqs(:)) + 1])
ylabel({'Virtual Range (km)'})
xlabel('Tx Freq (MHz)')
title(sprintf('Tx Loc: %1.1f N %1.1f E %1.1f km\n Mode: %i', ...
    rays(1).txloc(1), rays(1).txloc(2), rays(1).txloc(3), rays(1).OX_mode))
set(gca, 'YDir','reverse')
clim(climit)
hc = colorbar;
hc.Label.String = 'Received elevation (°)';
% hc.Label.String = 'Received azimuth (°)';
% hc.Label.String = 'Total absorption (dB)';
%legend({'X', 'No B', 'O'})
% title(sprintf('foF2: %1.1f MHz, hmF2: %1.1f km, S/C alt: %1.1f km', D.fof2(timeidx), D.hmf2(timeidx), alt_iss))

%set(gca,'CLim',colorlimit)

hold off

grid on
grid minor
