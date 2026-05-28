function range = raytrace_1d(alts_in, freqs, Ne, B, Inc, OX_mode, txalt)
%% raytrace_1d
% % 1D 'raytracer', returns reflection heights for specified frequencies
% % and ionosphere
% % Example
% 
% 
% alts = 90:2:460;
% txalt = 580
% lat = 20;
% lon = 0;
% year = 2020;
% [~, ~, ~, I, B] = igrfmagm(alts, lat * ones(size(alts)), ...
%     lon * ones(size(alts)), year * ones(size(alts)));
% Inc = I';  % inclination
% B = B' / 1E9;  % field strength
% Np = 1E12; hp = 300; H = 100;
% Ne = calc_chapman_prof(alts, Np, hp, H);
% ranges = raytrace_1d(alts, freqs, Ne, B, Inc, 1);
% clf
% subplot(1, 2, 1)
% hold on 
% plot(freqs, ranges)
% xlabel('Wave Freq (MHz)')
% legend({'O-mode'})
% grid on
% grid minor
% hold off
% subplot(1, 2, 2)
% plot(elec2freq(Ne) /1E6, alts)
% xlabel('Plasma Freq (MHz)')
% grid on
% grid minor


%% check if we're topside sounding 
hmf2 = alts_in(Ne == max(Ne));
if txalt > hmf2
    ts_sound = true;
end


%% Reorder profile and B/I to have transmitter at 0 alt
alts = -(alts_in - txalt);
idx = alts >= 0;
Ne = Ne(idx);
alts = alts(idx);
B = B(idx);
Inc = Inc(idx);

if ts_sound
    B = flip(B);
    Inc = flip(Inc);
    alts = flip(alts);
    Ne = flip(Ne);
end


%% % interpolate for the 'raytracer'
alt_stepsize = 0.001;
alti = min(alts):alt_stepsize:max(alts);
Ne_i = interp1(alts, Ne, alti);
I_i = interp1(alts, Inc, alti);
B_i = interp1(alts, B, alti);


%% Loop over freqs
range = zeros(size(freqs));
time = 0;
theta = deg2rad(I_i - 90);

for f = 1:length(freqs)
    freq = freqs(f);
    %% Calculate refractive index
    n = appleton_hartree_fast(theta, Ne_i, B_i, freq * 1E6 * ones(size(I_i)), OX_mode);
    dist_to_iono = min(alti);
    
    %% Calculate distance
    idx = find(real(n) == 0, 1, 'first');
    alt_stepsize = mean(unique(diff(alti)));  % km
    if isempty(idx)
        range(f) = NaN;
    else
        range(f) = sum(alt_stepsize ./ (n(1:idx - 1))) + dist_to_iono;
    end
end

















