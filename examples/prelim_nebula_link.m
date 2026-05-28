clear
close all

%% inputs
Pt = 50; % watts
N0_thermal = -174;
N0_above_thermal = 40;

t_int = 0.1;
polarization_loss = 3;
ionospheric_loss = 3;

freq_lims = [2e6 40e6];
range_lims = [200e3, 3000e3];

%% calcs
coeffs = polyfit(range_lims, freq_lims, 1);

ranges = linspace(range_lims(1), range_lims(2), 100);

% read ant file
ant_fname = 'ex_ant_gain.txt';
fid = fopen(ant_fname, 'r');
cnt = 1;
while ~feof(fid)

    raw = fgetl(fid);
    str_cell = strsplit(raw, ' ');
    ant_freq(cnt) = str2double(str_cell{2});
    ant_gain(cnt) = str2double(str_cell{3});
    cnt = cnt + 1;

end
fclose(fid);

for k = 1:length(ranges)

    r = ranges(k);
    
    f = r * coeffs(1) + coeffs(2);
    meas_freqs(k) = f;
    [~, gain_idx] = min(abs(f / 1e6 - ant_freq));
    G = ant_gain(gain_idx);

    N0 = N0_thermal + N0_above_thermal;
    eirp = 10*log10(Pt) + G; 
    
    % assume "infinite plane" for ionosphere
    fspl = 20*log10(r * f * 4*pi / 3e8);
    p_rx = eirp - fspl + G - polarization_loss - ionospheric_loss;
    
    c_n0 = p_rx - N0;
    coh_gain = 10*log10(t_int);
    snr(k) = c_n0 + coh_gain;

end

%% plotting
figure; 
plot(ranges/1000, snr, 'linewidth', 2);
grid;
yline(12, '--', 'Threshold', 'Color', 'r')
xlabel('Path Range (km)');
ylabel('Pulse SNR (dB)')
title('Pulse SNR vs Path Range')
%{
figure; 
plot(ranges / 1000, meas_freqs / 1e6, 'linewidth', 2)
grid;
xlabel('Path Range (km)')
ylabel('Measurement Frequency (MHz)')
title('Measurement Frequency vs Path Range')
%}
