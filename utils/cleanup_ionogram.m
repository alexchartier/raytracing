function [ionogram, idx] = cleanup_ionogram(homed_rays)
%% cleanup_ionogram
% Function to retrieve the 'good' part of an ionogram
%
% % Example
%
% mod_fn_fmt = '/Users/chartat1/data/sami3/2017_tid/sami_mat/{YYYY-mm-dd_HHMM}.mat';
% ray_fn_fmt = '/Users/chartat1/data/sami3/2017_tid/recon_rays/sim_twosat_rays_%s_%i.mat';
% 
% time = datenum(2017, 1, 10, 18, 0, 0);
% sami = loadstruct(filename(mod_fn_fmt, time));
% homed_rays = loadstruct(sprintf(ray_fn_fmt, 'vert', 2));
%
% [ionogram, idx] = cleanup_ionogram(homed_rays);
%
% % plotting to confirm we have a 'good' ray
% subplot(1, 2, 1)
% for i = 1:length(homed_rays{idx})
%     if homed_rays{idx}(i).initial_bearing == ionogram.az(idx) & ...
%             homed_rays{idx}(i).initial_elev == ionogram.el(idx)
%
%         nemax = max(homed_rays{idx}(i).electron_density);
%         nemaxid = homed_rays{idx}(i).electron_density == nemax;
%         lat = homed_rays{idx}(i).lat(nemaxid);
%         lon = homed_rays{idx}(i).lon(nemaxid);
%         height = homed_rays{idx}(i).height(nemaxid);
%     end
% end
% lon(lon < 0) = lon(lon < 0) + 360;
%
% ne = zeros(size(sami.alt));
% for i = 1:length(ne)
%     ne(i) = interp2(sami.lat, sami.lon, squeeze(sami.dene(i, :, :))', lat, lon);
% end
% hold on
% plot(elec2freq(ne) / 1E3, sami.alt, '-k', 'LineWidth', 3)
% plot(elec2freq(nemax) / 1E3, height, 'r.', 'MarkerSize', 40)
% hold off
% xlabel('Freq (MHz)')
% ylabel('Alt (km)')
% grid on; grid minor
% 
% subplot(1, 2, 2)
% plot(ionogram.freqs, ionogram.vht)
% xlabel('Freq (MHz)')
% ylabel('Virtual height (km)')
% grid on; grid minor

%% Check if we're vertical or oblique
if sum(homed_rays{1}(1).txloc == homed_rays{1}(1).rxloc) == 3
    vert = true;
else
    vert = false;
end
branched = false;

%% Figure out the maximum number of rays 
maxrays = 0;
freqs = [];
for ri = 1:length(homed_rays)
    freqs = [freqs, homed_rays{ri}(1).frequency];
    if length(homed_rays{ri}) > maxrays
        maxrays = length(homed_rays{ri});
    end
end

%% Go through and get all the az/el at each freq
az = zeros(length(freqs), maxrays) * NaN;
el = zeros(length(freqs), maxrays) * NaN;
rg = zeros(length(freqs), maxrays) * NaN;
lat_ref = zeros(length(freqs), maxrays) * NaN;
lon_ref = zeros(length(freqs), maxrays) * NaN;
for f = 1:length(freqs)
    for r = 1:length(homed_rays{f})
        az(f, r) = homed_rays{f}(r).initial_bearing;
        el(f, r) = homed_rays{f}(r).initial_elev;
        rg(f, r) = homed_rays{f}(r).group_range_to_rx;

        hts = homed_rays{f}(r).height;
        hti = hts == min(hts);
        lat_ref(f, r) = homed_rays{f}(r).lat(hti);
        lon_ref(f, r) = homed_rays{f}(r).lon(hti);
    end
end

%% For each frequency, select only the lowest-elevation ones
az1 = zeros(size(freqs));
el1 = zeros(size(freqs));
rg1 = zeros(size(freqs));
lat_r1 = zeros(size(freqs));
lon_r1 = zeros(size(freqs));
for f = 1:length(freqs)
    idx = el(f, :) == nanmin(el(f, :));
    az1(f) = az(f, idx);
    el1(f) = el(f, idx);
    rg1(f) = rg(f, idx);
    lat_r1(f) = lat_ref(f, idx);
    lon_r1(f) = lon_ref(f, idx);
end

%% Discard the last part of the ionogram if there's a break in vector angle
theta = vec_angle(az1(1:end-1), el1(1:end-1), az1(2:end), el1(2:end));
idx = find(abs(theta) > 15);
if idx
    freqs = freqs(1:idx);
    rg1 = rg1(1:idx);
    az1 = az1(1:idx);
    el1 = el1(1:idx);
    lat_r1 = lat_r1(1:idx);
    lon_r1 = lon_r1(1:idx);

    if vert
        branched = true;
    end
end

% %% Discard the top end of the ionogram if it looks bad
% id2 = find(diff(rg1) < -50);
% 
% if id2
%     idx = min(id2);
%     freqs = freqs(1:idx);
%     rg1 = rg1(1:idx);
%     az1 = az1(1:idx);
%     el1 = el1(1:idx);
%     lat_r1 = lat_r1(1:idx);
%     lon_r1 = lon_r1(1:idx);
%     branched = true;
% end


%% Get params
OX_mode = homed_rays{1}(1).OX_mode;
N_local = homed_rays{1}(1).electron_density(1);

ionogram.freqs = freqs;
ionogram.vht = rg1 / 2;
ionogram.az = az1;
ionogram.el = el1;
ionogram.mode = OX_mode;
ionogram.N_local = N_local;
ionogram.txloc = homed_rays{1}(1).txloc; 
ionogram.lat_ref = lat_r1;
ionogram.lon_ref= lon_r1;
ionogram.branched = branched;


end

%% vec_angle
function angle = vec_angle(az1, el1, az2, el2)
% Convert az-el vectors to Cartesian coordinates
cart1 = sphcart([ones(size(az1)); deg2rad(el1); deg2rad(az1)]');
cart2 = sphcart([ones(size(az2)); deg2rad(el2); deg2rad(az2)]');

% Calculate the angle using the dot product
angle = rad2deg(acos(dot(cart1',  cart2')));
end








