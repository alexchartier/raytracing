%% twosat_synth.m
% Simulate two satellite topside sounder mission data

clear

%% inputs
out_fn = 'sim_twosat_rays.mat';
time = datenum(2019, 9, 21, 0, 0, 0);
freqs = 3:0.1:15; % 2:0.2:20; % 20;

alts = 150:2:800;
sat_alt = 600;

satlats = [47, 52];
satlons = [-75];
satlons(satlons > 180) = satlons(satlons > 180) - 360;
[satlatarr, satlonarr] = meshgrid(satlats, satlons);
satlocs = [satlatarr(:), satlonarr(:), ones(numel(satlatarr), 1) * sat_alt];

% grid
lats = 30:2:65;
lons = -80:2:-70;
lons(lons > 180) = lons(lons > 180) - 360;

elvarr = [-90:5:-30];
azarr = -180:5:180;
R12 = 100;

OX_mode = 0;

kp = 3;
maxdist = 1E5;  % meters from homing
blind_range = 30E3; % two-way (e.g. there-and-back range)
tol = [1e-7 0.01 25];       % ODE solver tolerance and min max stepsizes



%% Load ionosphere
[iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    Bx, By, Bz, geomag_grid_parms] = ...
    gen_iono_geomag_grids(alts, lats, lons, time, R12);

%% Loop over transmitters
clear homed_rays
fprintf('Started a new raytracing expt\n')
for s1 = 1%:size(satlocs, 1)
    txloc = satlocs(s1, :);
    %% loop over receivers
    for s2 = 2%1:size(satlocs, 1)
        rxloc = satlocs(s2, :);


        %% Get plasma freq and refractive index at the txloc
        iono_pf_grid = sqrt(80.6 * iono_en_grid ./ 1E6);
        plasmafreq = interp2(lats, lons, squeeze(iono_pf_grid(:, :, alts == sat_alt))', txloc(1), txloc(2));
        assert(length(plasmafreq) == 1, 'Come back and fix this - make a real interpolator')

        iono_pf_prof = nan(size(alts));
        for i = 1:length(alts)
            iono_pf_prof(i) = interp2(lats, lons, squeeze(iono_pf_grid(:, :, alts == alts(i)))', txloc(1), txloc(2));

        end
        fprintf('FoF2: %1.2f\n', max(iono_pf_prof))


        %% Raytrace each frequency that's above the local plasma frequency
        for f = 1:length(freqs)
            fprintf('Freq: %1.1f MHz\n', freqs(f))
            %% Local plasma > transmission freq - skip
            if plasmafreq > freqs(f)
                %disp('local plasma environment > transmit freq - skipping')
                continue
            end

            %% Refractive index calculation.
            % Estimate the refractive  index at ray origin. Note this is only the
            % formula for the no mag. field case. Recommend the use of the full
            % Appleton-Hartree equation to calculate the refractive index for the
            % O and X modes
            refractive_ind = sqrt(1 - plasmafreq .^ 2 ./ freqs(f) .^ 2);
            
            fprintf('Refractive Index: %1.3f \n', refractive_ind)


            %% Raytrace the problem

            homed_ray = ...
                raytrace_itsi(freqs(f), OX_mode, txloc, rxloc, iono_en_grid, iono_en_grid_5, ...
                collision_freq, iono_grid_parms, Bx, By, Bz, geomag_grid_parms, ...
                elvarr, azarr, maxdist, tol, refractive_ind, blind_range);

            if isstruct(homed_ray)
                homed_rays(s1, s2, f) = homed_ray;
            end

        end
    end
end
clear raytrace_3d

% save
homed_rays(1).iono_en_grid = iono_en_grid;
savestruct(out_fn, homed_rays)
fprintf("saved homed_rays to %s\n", out_fn)

%% load

homed_rays = loadstruct(out_fn);
iono_en_grid = homed_rays(1).iono_en_grid;



%% Plotting ionospheric density slice
clf


% convert the coodinate frame to curved Earth geometry
[lon3, lat3, alt3] = meshgrid(lons(2), lats, alts);

CartV = sphcart([(6371 + alt3(:)) * 1E3, deg2rad(lat3(:)), deg2rad(lon3(:))]);

V = iono_en_grid(:, 2, :);
X = reshape(CartV(:, 1), size(V));
Y = reshape(CartV(:, 2), size(V));
Z = reshape(CartV(:, 3), size(V));

hold on
satglobe4e
iono_panel = surf(squeeze(X), squeeze(Y), squeeze(Z), 'FaceColor', 'none', 'EdgeColor', 'none');
alpha   = 0.7; % globe transparency level, 1 = opaque, through 0 = invisible
set(iono_panel, 'FaceColor', 'texturemap', 'CData', squeeze(V), 'FaceAlpha', alpha, 'EdgeColor', 'none');

% plot the rays
rays = homed_rays(:);
Re = 6371E3;
for r = 1:length(rays)
    if ~isempty(rays(r).initial_elev)
        hidx = rays(r).height <= rxloc(3);
        sph = [rays(r).height(hidx) * 1E3 + Re; deg2rad(rays(r).lat(hidx)); deg2rad(rays(r).lon(hidx))];
        cart = sphcart(sph');
        plot3(cart(:, 1), cart(:, 2), cart(:, 3), 'w')
    end
end

for l = 1:size(satlocs, 1)
    loc = satlocs(l, :);
    cart = sphcart([loc(3) * 1E3 + Re, deg2rad(loc(1)), deg2rad(loc(2))]);
    plot3(cart(1), cart(2), cart(3), 'ro', 'markersize', 10, 'markerfacecolor', 'r')
end

hold off

%%

set(gca,'color','black')
export_fig('~/Downloads/twosat_sim.tif');


%% Plotting 2 - ionogram plot

figure
colormap('cool')

xlimit = [2, 12];
ylimit = [0, 3000];
colorlimit = [-130 -80];
ct = 0;
hs = [];
for s1 = 1:size(homed_rays, 1)
    for s2 = s1:size(homed_rays, 2)
        ct = ct + 1;
        hs = [hs, subplot(size(homed_rays, 2), 1, ct)];

        rays = squeeze(homed_rays(s1, s2, :));
        hold on

        for r = 1:length(rays)
            if rays(r).home
                freq = rays(r).frequency;
                rg = rays(r).group_range_to_rx;
                loss = rays(r).absorption(end) + fspl(rg * 1E3, 3E8 / (freq * 1E6));
                scatter(freq, rg, 150, -loss, 'filled', 'o')
            end
        end

        ylim(ylimit)
        xlim(xlimit)
        ylabel({'Virtual Range (km)'})
        if ct == size(homed_rays, 2)
            xlabel('Tx Freq (MHz)')
        else
            set(gca, 'XTickLabels', [])
        end
        set(gca, 'YDir','reverse')
        %legend({'X', 'No B', 'O'})
        % title(sprintf('foF2: %1.1f MHz, hmF2: %1.1f km, S/C alt: %1.1f km', D.fof2(timeidx), D.hmf2(timeidx), alt_iss))

            set(gca,'CLim',colorlimit)

        hold off

        grid on
        grid minor
        
    end
end


for ct = 1:length(hs)
        hp = get(hs(ct),'Position');
        set(hs(ct), 'Position', [hp(1), hp(2), hp(3) * 0.8, hp(4)])
        
end
hp3 = get(hs(end),'Position');
 
h = colorbar('Position', [hp3(1)+hp3(3)+0.02  hp3(2)+0.05  0.02  hp3(2)+hp3(3)]);
h.Label.String = 'Signal Loss (dB)';
    set(gcf,'Color','w')
    set(gca, 'FontSize', 20, 'CLim',colorlimit)
        %set(h, 'FontSize', 20)
        
    set(findall(gcf,'-property','FontName'),'fontsize',18) 
        











