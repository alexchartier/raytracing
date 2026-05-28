function rg = vis_synth_chap(h_peak, N_peak, N_local, freqs, OX_mode, alts, Inc, B)

%% Simulate vertical ionograms from Chapman ionosphere

%% increment H to get a match to local Ne
H = 1;
Ne_out = 0;
Ne = calc_chapman_prof(alts, N_peak, h_peak, H, c);
while N_local > Ne(1)
    Ne = calc_chapman_prof(alts, N_peak, h_peak, H);
    H = H + 1;
end
rg = raytrace_1d(alts, freqs, Ne, B, Inc, OX_mode);

disp(H)

%% return 


