%% Simulate vertical ionograms from IRI
% % Group velocity derivation
% % n = c / vp; 
% % 
% % vp * vg = c^2; 
% % vp = c^2 / vg;
% % 
% % n = vg / c; 
% % vg = n * c;
% %
% % Example:

% alts = 200:2:580;
%%
UT = [2018 1 15 0 0];         % UT - year, month, day, hour, minute
R12 = 1;                    % yearly smoothed sunspot number\

foE = 1;
foF2 = 6;
hmF2 = 240;
freqs = 3:0.1:12; 

lat = 38;              % latitude of the start point of ray (Darwin)
lon = 0;    
alts = 0:2:600;
alti = 200:0.001:580;

ts_sound = true;


%% Run IRI and IGRF
iri_options = [];
iri_options.foE = foE;  % 2 - 6 good 
iri_options.foF2 = foF2;  % 6 - 10
iri_options.hmF2 = hmF2;  % 250:25:350

% run IRI
[iono, iono_extra] = iri2020(lat, lon, R12, UT, alts(1), ...
    alts(2) - alts(1), length(alts), iri_options);
Ne = iono(1, :);


% run IGRF
[~, ~, ~, I, B] = igrfmagm(alts, lat * ones(size(alts)), ...
    lon * ones(size(alts)), UT(1) * ones(size(alts))); 
I = I';  % inclination
B = B' / 1E9;  % field strength

Ne_i = interp1(alts, Ne, alti);
I_i = interp1(alts, I, alti);
B_i = interp1(alts, B, alti);


%% Loop over freqs
O_mode = zeros(size(freqs));
X_mode = zeros(size(freqs));
No_B = zeros(size(freqs));
for f = 1:length(freqs)
    freq = freqs(f);
    %% Calculate refractive index
    theta = deg2rad(I_i - 90);
    [n_O, n_X, n] = appleton_hartree(theta, Ne_i, B_i, freq * 1E6 * ones(size(I_i)));

    if ts_sound
        n_O = n_O(end:-1:1);
        n_X = n_X(end:-1:1);
        n = n(end:-1:1);
        dist_to_iono = 0;
    else
        dist_to_iono = min(alti);
    end

    
    %% Calculate distance
    O_idx = find(real(n_O) == 0, 1, 'first');
    X_idx = find(real(n_X) == 0, 1, 'first');
    No_B_idx = find(real(n) == 0, 1, 'first');
    alt_stepsize = mean(unique(diff(alti)));  % km
    O_mode(f) = sum(alt_stepsize ./ (n_O(1:O_idx - 1))) + dist_to_iono;
    X_mode(f) = sum(alt_stepsize ./ (n_X(1:X_idx - 1))) + dist_to_iono;
    No_B(f) = sum(alt_stepsize ./ (n(1:No_B_idx - 1))) + dist_to_iono;

end
O_mode(O_mode == min(alts)) = NaN;
X_mode(X_mode == min(alts)) = NaN;


%% 
clf
subplot(1, 2, 1)
hold on 
plot(freqs, O_mode)
plot(freqs, X_mode)
grid on
grid minor
hold off
subplot(1, 2, 2)
plot(elec2freq(Ne_i) /1E6, alti)
grid on
grid minor

% %%
% hold on
% % plot(elec2freq(Ne) / 1E6, alts)
% plot(real(n_O), alts)
% plot(real(n_X), alts)
% plot(real(n), alts)
% legend({'O', 'X', 'No B'})
% 
% hold off
% 



