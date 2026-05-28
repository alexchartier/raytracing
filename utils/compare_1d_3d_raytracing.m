%% compare_1d_3d_raytracing.m
% check that the 1D 'raytracer' is doing something similar to the 3D
% raytracer
clear
%% inputs
time = datenum(2017, 1, 10, 18, 0, 0);
mod_fn_fmt = '/Users/chartat1/data/sami3/2017_tid/sami_mat/{YYYY-mm-dd_HHMM}.mat';
ray_fn_fmt = '/Users/chartat1/data/sami3/2017_tid/recon_rays/sim_twosat_rays_%s_%i.mat';
freqs = 2:0.1:7;

%% load
sami = loadstruct(filename(mod_fn_fmt, time));
homed_rays = loadstruct(sprintf(ray_fn_fmt, 'vert', 2));
[ionogram, idx] = cleanup_ionogram(homed_rays);


%% figure out where to look
for i = 1:length(homed_rays{idx})
    if homed_rays{idx}(i).initial_bearing == ionogram.az(idx) & ...
            homed_rays{idx}(i).initial_elev == ionogram.el(idx)

        nemax = max(homed_rays{idx}(i).electron_density);
        nemaxid = homed_rays{idx}(i).electron_density == nemax;
        lat = homed_rays{idx}(i).lat(nemaxid);
        lon = homed_rays{idx}(i).lon(nemaxid);
        height = homed_rays{idx}(i).height(nemaxid);
    end
end
lon(lon < 0) = lon(lon < 0) + 360;

lat = 38; lon = 282.5;

Ne = zeros(size(sami.alt));
for i = 1:length(Ne)
    Ne(i) = interp2(sami.lat, sami.lon, squeeze(sami.dene(i, :, :))', lat, lon);
end
Ne = Ne * 1E6;
alts = sami.alt;
txalt = homed_rays{1}(1).txloc(3);
OX_mode = homed_rays{1}(1).OX_mode;


%% calculate a 1D ionogram on the same density grid
[~, ~, ~, Inc, B] = igrfmagm(alts, lat * ones(size(alts)), ...
    lon * ones(size(alts)), str2num(datestr(time, 'YYYY')) * ones(size(alts)));
Inc = Inc';  % inclination
B = B' / 1E9;  % field strength
vht = raytrace_1d(alts, freqs, Ne, B, Inc, OX_mode, txalt);


%% plotting
subplot(1, 2, 1)

hold on
plot(elec2freq(Ne) / 1E3, sami.alt, '-k', 'LineWidth', 3)
plot(elec2freq(nemax) / 1E3, height, 'r.', 'MarkerSize', 40)
hold off
xlabel('Freq (MHz)')
ylabel('Alt (km)')
grid on; grid minor

subplot(1, 2, 2)
hold on
plot(ionogram.freqs, ionogram.vht, '-k', 'LineWidth', 3)
plot(freqs, vht, '-r', 'LineWidth', 3)
hold off
xlabel('Freq (MHz)')
ylabel('Virtual height (km)')
legend({'3D', '1D'})
grid on; grid minor







