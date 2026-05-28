%% plot the simulated topside ionograms
in_fn_fmt = './../data/iri_osse_el80-0_az0-360/sim_twosat_rays_O_1.mat';
out_fn_fmt = 'data/plots/sim_twosat_rays_%i.png';

modes = {'O', 'X'};
for i = 1 % :20
    clf
    hold on
    plot(NaN, NaN, 'g.', 'MarkerSize', 30)
    plot(NaN, NaN, 'r.', 'MarkerSize', 30)
    for m = modes
        switch m{1}
            case 'O'
                plt = 'g.';
            case 'X'
                plt = 'r.';
        end
        homed_rays = loadstruct(sprintf(in_fn_fmt, m{1}, i));

        for r = 1:length(homed_rays)
            for l = 1:length(homed_rays{r})
                plot(homed_rays{r}(l).frequency, ...
                    homed_rays{r}(l).group_range_to_rx, plt, 'MarkerSize', 10)
            end
        end
    end
    title(sprintf('%i km alt, Sat A: %1.1f N, Sat B: %1.1f N', ...
        round(homed_rays{r}(l).height(1)), homed_rays{r}(l).lat(1), homed_rays{r}(l).lat(1)+ 5))
    set(gca, 'YDir', 'reverse', 'FontSize', 20)
    xlabel('Freq (MHz)')
    ylabel('Group Range (km)')
    ylim([500, 2000])
    xlim([5, 12])
    legend(modes)
    grid on
    grid minor
    %export_fig(sprintf(out_fn_fmt, i))
end

%% Plotting ionospheric density slice
homed_rays = loadstruct(sprintf(in_fn_fmt, 'O', 10));

iono_en_grid = homed_rays{1}(1).iono_en_grid;
lats = homed_rays{1}(1).iono_lat;
lons = homed_rays{1}(1).iono_lon;
alts = homed_rays{1}(1).iono_alt;
rxloc = homed_rays{1}(1).rxloc;

close
% convert the coodinate frame to curved Earth geometry
[lon3, lat3, alt3] = meshgrid(lons, lats, alts);

CartV = sphcart([(6371 + alt3(:)) * 1E3, deg2rad(lat3(:)), deg2rad(lon3(:))]);

V = iono_en_grid(:, :, :);
X = reshape(CartV(:, 1), size(V));
Y = reshape(CartV(:, 2), size(V));
Z = reshape(CartV(:, 3), size(V));

hold on
satglobe4e
%earth_example
iono_panel = surf(squeeze(X), squeeze(Y), squeeze(Z), 'FaceColor', 'none', 'EdgeColor', 'none');
alpha   = 0.7; % globe transparency level, 1 = opaque, through 0 = invisible
set(iono_panel, 'FaceColor', 'texturemap', 'CData', squeeze(V), 'FaceAlpha', alpha, 'EdgeColor', 'none');

% plot the rays
% rays = homed_rays(:);
Re = 6371E3;
for ri = 1:length(homed_rays)
    rays = homed_rays{ri};
    for r = 1:length(rays)
        hidx = rays(r).height <= rxloc(3);
        sph = [rays(r).height(hidx) * 1E3 + Re; deg2rad(rays(r).lat(hidx)); deg2rad(rays(r).lon(hidx))];
        cart = sphcart(sph');
        plot3(cart(:, 1), cart(:, 2), cart(:, 3), 'w')
    end
end

satlocs = zeros(2, 3);
satlocs(1, :) = [homed_rays{1}(1).lat(1), homed_rays{1}(1).lon(1), homed_rays{1}(1).height(1)];
satlocs(2, :) = homed_rays{1}(1).rxloc;
for l = 1:size(satlocs, 1)
    loc = satlocs(l, :);
    cart = sphcart([loc(3) * 1E3 + Re, deg2rad(loc(1)), deg2rad(loc(2))]);
    plot3(cart(1), cart(2), cart(3), 'ro', 'markersize', 10, 'markerfacecolor', 'r')
end

      R = 140*6400E3; 

            set(gca,'CameraViewAngle',.85,'CameraPosition',[1,0,0]*R,'CameraUpVector',[0,0,1]);

hold off

%%
set(gca,'color','black')
%export_fig('~/Downloads/twosat_sim.tif');
