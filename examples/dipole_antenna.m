%% simulate dipole antenna performance

%% Work out what the optimal impedance is, looking at the low end

%% inputs/outputs
len = 12;
wdt = 0.20;

Zs = 50:100:4000;
freqs = [2, 6, 12, 20] * 1E6; 

out_fn_pattern = '~/Documents/Papers/2024_nebula/plots/antenna/%1.1f_m_%1.1f_MHz_pattern.png';
out_fn_impedance = '~/Documents/Papers/2024_nebula/plots/antenna/%1.1f_m_%1.1f_MHz_impedance.png';
out_fn_gain_vs_freq = '~/Documents/Papers/2024_nebula/plots/antenna/%1.1f_m_gain_vs_freq.png';


%% Set up antenna
radius = strip2cylinder(wdt);
d = dipoleCylindrical ('length', len, 'radius', radius);

%% loop over freqs
for f = 1:length(freqs)
    %%
    close all    
    freq = freqs(f);

    titlestr = sprintf('%i m tip-to-tip dipole, Freq: %i MHz', len, freq/1E6);

    radius = strip2cylinder(wdt);
    d.TiltAxis = 'X';
    d.Tilt = 90;
    realized_gains = zeros(size(Zs));

    for zi = 1:length(Zs)
        Z = Zs(zi);
        s = sparameters(d,freq, Z);

        S11 = squeeze(s.Parameters);

        realized_gains(zi) = (1 - abs(S11).^2) .* max(10.^(patternAzimuth(d, freq)/10));
    end

    realized_gain_db = real(10 * log10(realized_gains));
    max_Z = Zs(realized_gain_db == max(realized_gain_db));
    
    % Antenna pattern plot
    figure('Position',[600 100 1200 800])
    pattern(d, freq)
    title(titlestr)
    set(gca, 'FontSize', 24, 'FontName', 'Futura')
    saveas(gcf, sprintf(out_fn_pattern, len, freq/1E6))

    % Impedance plot
    figure('Position',[600 100 1200 800])
    plot(Zs, realized_gain_db, '-x')
    title(sprintf('%1.1f m dipole, %1.1f MHZ,  Max gain: %1.1f dB @ %i ohms', ...
        len, freq/1E6, max(realized_gain_db), max_Z))

    xlabel('Input Impedance (Ohms)')
    ylabel('Free Space Realized Gain (dBi)')
    
    set(gca, 'FontSize', 24, 'FontName', 'Futura')
    grid on
    grid minor

    saveas(gcf, sprintf(out_fn_impedance, len, freq/1E6))


end

%% figure out the maximum gain if using antenna matched to the lowest frequency
Z = 1650; 
realized_gains = zeros(size(freqs));
for f = 1:length(freqs)
    freq = freqs(f);
    s = sparameters(d, freq, Z);
    S11 = squeeze(s.Parameters);
    realized_gains(f) = (1 - abs(S11).^2) .* mean(patternAzimuth(d, freq));

end
plot(freqs/1E6, 10 * log10(realized_gains), '-x', 'LineWidth', 3);
xlabel('Freq (MHz)')
ylabel(sprintf('Realized Gain w. %1.1f-m dipole @ %i ohms', len, Z))
grid on
grid minor
set(gca, 'FontSize', 24, 'FontName', 'Futura')
% saveas(gcf, sprintf(out_fn_gain_vs_freq, len))


%%
Z = 1650; 

freq = 3E6;
lengths = 1:10;
radius = strip2cylinder(wdt);
realized_gains = zeros(size(lengths));

for l = 1:length(lengths)
    d = dipoleCylindrical ('length', lengths(l), 'radius', radius);
    s = sparameters(d, freq, Z);
    S11 = squeeze(s.Parameters);
    realized_gains(l) = (1 - abs(S11).^2) .* mean(patternAzimuth(d, freq));
end

plot(lengths, 10 * log10(realized_gains), '-x')
xlabel('Dipole Length (m)')
ylabel('Gain @ 3 MHz (dB)')
grid on
grid minor


%% Plot the Smith charts
%Z = 50;
% figure
% s = sparameters(d,freq, Z);
% s_off = sparameters(d_off,freq, Z);
% clf
% subplot(1, 2, 1)
% hg = smithplot(s,1,1,'GridType','ZY');
% hg.LineStyle = '--';
% hg.TitleBottom = sprintf('%i-m dipole, %i ohm ',len, Z);
% 
% 
% subplot(1, 2, 2)
% hg = smithplot(s_off,1,1,'GridType','ZY');
% hg.LineStyle = '--';
% hg.TitleBottom = sprintf('%i-m dipole, feed offset by %i-m, %i ohm ', len, offset, Z);
% 
% 
