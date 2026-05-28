function prof = gen_iri_profile(time, NmF2, hmF2, N_sc, h_sc, lat, lon, alts, R12)
%% gen_iri_profile.m
% 
% time = datetime(2017, 1, 1);
% NmF2 = 4E5; 
% hmF2 = 300;
% N_sc = 1E4; 
% h_sc = 580;
% 
% R12 = 30;
% lat = 38;
% lon = -77.5;
% alts = 90:2:600;
% prof = gen_iri_profile(time, NmF2, hmF2, N_sc, h_sc, lat, lon, alts, R12)
% 
% hold on
% plot(prof, alts)
% plot(N_sc, h_sc, 'xr', 'MarkerSize', 10)
% hold off

%%
UT = [year(time), month(time), day(time), hour(time), minute(time)];
iri_options.hmF2  = hmF2;
iri_options.foF2  = elec2freq(NmF2) / 1E3;
[iono, iono_extra] = iri2016(lat, lon, R12, UT, ...
            min(alts), alts(2) - alts(1), length(alts), iri_options);
prof = iono(1, :) / 1E6;%   + N;
N_sc_mod = interp1(alts, prof, h_sc);
scalefac = N_sc / N_sc_mod;
tsi = alts > hmF2; 
scaling = ((alts(tsi) - hmF2).^2 ./ (h_sc - hmF2).^2);
prof(tsi) = prof(tsi) .* (1 - scaling) + (prof(tsi) .* scalefac) .* scaling;

% function read_ig_rz(ig_rz_path)
% 
% ig_rz_path = '~/pharlap/dat/iri2016/ig_rz.dat';
% txt = asciiread(ig_rz_path);
% hdr = txt(1:4, :);
% txt = txt(5:end, :);
% 
% txtf = [];
% for i = 1:size(txt, 1)
%     txtf = [txtf, txt(i, :)];
% end
% rz = str2num(txtf);
% timestr = str2num(hdr(3, :));
% stime = datetime(timestr(2), timestr(1), 1);
% etime = datetime(timestr(4), timestr(3), 1);
% months = stime:calmonths(1):etime;

% end