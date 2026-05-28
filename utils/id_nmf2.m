function nmf2 = id_nmf2(fmax, OX_mode, alts, Inc, B)
%% get the nmf2 from the max freq of the ionogram


nmf2 = freq2elec(fmax) * 1E12; 

h_peak = max(alts);  % peak height (distance below s/c)
H_dummy = 100;
Ne = calc_chapman_prof(alts, nmf2, h_peak, H_dummy);
range = raytrace_1d(alts, fmax, Ne, B, Inc, OX_mode);

%% increment through to get correct Np
if isnan(range)
    while isnan(range)
        Ne = calc_chapman_prof(alts, nmf2, h_peak, H_dummy);
        range = raytrace_1d(alts, fmax, Ne, B, Inc, OX_mode);
        nmf2 = nmf2 + 1E9;
        %disp(nmf2)
    end
else
    while ~isnan(range)
        Ne = calc_chapman_prof(alts, nmf2, h_peak, H_dummy);
        range = raytrace_1d(alts, fmax, Ne, B, Inc, OX_mode);
        nmf2 = nmf2 - 1E9;
        %disp(nmf2)
    end
end

fprintf('NmF2: %1.1e\n', nmf2)
