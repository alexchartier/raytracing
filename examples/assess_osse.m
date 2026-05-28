%% assess_osse.m
% title_str = '750-km spacing case, low & high latitude rays';
% fn2 = '~/data/nebula/iri_osse_el80-0_az0-360_2sats_hilat/homings/sat1_sat2/sim_twosat_rays_O_%i.mat';
% fn = '~/data/nebula/iri_osse_el80-0_az0-360_2sats_lowlat/homings/sat1_sat2/sim_twosat_rays_O_%i.mat';

title_str = '10-km spacing case, low & high latitude rays';
fn = '~/data/nebula/iri_osse_el80-0_az0-360_2sats_nvis_lowlat/homings/sat1_sat2/sim_twosat_rays_O_%i.mat';
fn2 = '~/data/nebula/iri_osse_el80-0_az0-360_2sats_nvis_hilat/homings/sat1_sat2/sim_twosat_rays_O_%i.mat';

%% load
D = loadstruct(sprintf(fn, 1));
D2 = loadstruct(sprintf(fn2, 285));
D = cat(2, D, D2);

%% Pull az/el info out
freqs = [];
geo_dist = [];
azimuths = [];
elevs = [];
OX_mode = [];
for ct = 1:length(D)
    txloc = [D{ct}(1).lat(1), D{ct}(1).lon(1), D{ct}(1).height(1)];
    sc_az = azimuth(D{ct}(1).lat(1), D{ct}(1).lon(1), D{ct}(1).rxloc(1), D{ct}(1).rxloc(2));

    freqs = [freqs; vertcat(D{ct}.frequency)];
    OX_mode = [OX_mode; vertcat(D{ct}.OX_mode)];
    geo_dist = [geo_dist; vertcat(D{ct}.geometric_dist_to_rx)];
    azimuths = [azimuths; vertcat(D{ct}.initial_bearing) - sc_az];

    elevs = [elevs; vertcat(D{ct}.initial_elev)];
end

%% histos
close
figure('Position', [10, 10, 1200, 600])
B = uicontrol(gcf, 'Style', 'Text', 'String', title_str, ...
    'Position', [200, 580, 800, 30], 'FontSize', 24, 'FontName', 'Futura')

% freq/range
subplot(1, 2, 1)
histogram2(geo_dist, freqs, [0:50:2000], 2:20, 'FaceColor', 'flat')
ylabel('Freq (MHz)')
xlabel('Geometric range (km)')

% clim(gca, [0, 1800])
y = colorbar;
ylabel(y, '# rays', 'rotation', 270, 'FontSize', 24)
view([0, 90])
grid on; grid minor

set(gca, 'FontSize', 24, 'FontName', 'Futura')

% Az/El hist
subplot(1, 2, 2)
azimuths(azimuths < -150) = azimuths(azimuths < -150) + 180;
azimuths(azimuths > 150) = azimuths(azimuths > 150) - 180;
histogram2(azimuths, elevs, -90:5:90, -90:5:0, 'FaceColor', 'flat')
ylabel('Initial Elevation (deg)')
xlabel('Initial Azimuth (deg clockwise of S/C X)')
y = colorbar;
ylabel(y, '# rays', 'rotation', 270, 'FontSize', 24)
grid on; grid minor
view([0, 90])
% clim(gca, [0, 10000])


set(gca, 'FontSize', 24, 'FontName', 'Futura')