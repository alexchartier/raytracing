%% constellation_synth.m
% Simulate topside sounder constellation mission data

clear

%% inputs
out_fn = 'sim_constell_rays.mat';
time = datenum(2019, 9, 21);
freqs = 5:20; % 2:0.2:20; % 20;

alts = 150:2:900;
sat_alt = 890;

satlats = [10:10:30];
satlons = [210:10:230];
satlons(satlons > 180) = satlons(satlons > 180) - 360;
[satlatarr, satlonarr] = meshgrid(satlats, satlons);
satlocs = [satlatarr(:), satlonarr(:), ones(numel(satlatarr), 1) * sat_alt];

% grid
lats = 8:32;
lons = 208:232;
lons(lons > 180) = lons(lons > 180) - 360;


elvarr = [-90:5:-20];
azarr = -180:5:180;
R12 = 100;

OX_mode = 0;

kp = 3;
maxdist = 1E5;  % meters from homing
blind_range = 30E3; % two-way (e.g. there-and-back range)
tol = [1e-7 0.01 25];       % ODE solver tolerance and min max stepsizes


%% Specify time in PHaRLaP format
year = str2double(datestr(time, 'yyyy'));
month = str2double(datestr(time, 'mm'));
day  = str2double(datestr(time, 'dd'));
hour = str2double(datestr(time, 'HH'));
minute = str2double(datestr(time, 'MM'));
UT = [year month day hour minute];


%% Load ionosphere

[iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
    Bx, By, Bz, geomag_grid_parms] = ...
    gen_iono_geomag_grids(alts, lats, lons, UT, R12);

%% Loop over transmitters
clear homed_rays
for s1 = 1:size(satlocs, 1)
    txloc = satlocs(s1, :);
    %% loop over receivers
    for s2 = 1:size(satlocs, 1)
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
            if plasmafreq > freqs(f)
                %disp('local plasma environment > transmit freq - skipping')
                continue
            end

            % Refractive index calculation.
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

%% Save/load

savestruct(out_fn, homed_rays)
homed_rays = loadstruct(out_fn);


%% Plotting
earth_example


%
% P.CameraTargetMode   = 'manual';
% P.CameraTarget       = [0,0,0];
% P.CameraPositionMode  = 'manual';
% P.CameraViewAngleMode = 'manual';
% set(gca,P);
% R = 140*6400E3;
% %set(gca,'CameraViewAngle',.85,'CameraPosition',[1,0,0]*R,'CameraUpVector',[0,0,1]);

rays = homed_rays(:);
hold on
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



%% Plotting 2 - ionogram plot


for s1 = 1:size(homed_rays, 1)
    for s2 = 1:size(homed_rays, 2)
        %%
                clf
        rays = squeeze(homed_rays(s1, s2, :));
        hold on

        for r = 1:length(rays)
            if rays(r).home
                freq = rays(r).frequency;
                rg = rays(r).group_range_to_rx;
                loss = rays(r).absorption(end) + fspl(rg* 1E3, 3E8/ (freq*1E6));
                scatter(freq, rg, 150, -loss, 'filled', 'o')
            end
        end

        ylim([0, 3000])
        xlim([0, max(freqs) + 2])
        ylabel('Group Range (km)')
        xlabel('Tx Freq (MHz)')
        set(gca, 'YDir','reverse')
        %legend({'X', 'No B', 'O'})
        % title(sprintf('foF2: %1.1f MHz, hmF2: %1.1f km, S/C alt: %1.1f km', D.fof2(timeidx), D.hmf2(timeidx), alt_iss))
        grid on
        grid minor
        hC = colorbar;
        hC.Label.String = 'Signal Loss (dBm)';
        set(gca, 'FontSize', 20, 'CLim', [-130 -100])
        set(hC, 'FontSize', 20)
        hold off
        pause(10)

    end
end













