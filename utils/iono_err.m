function err = iono_err(H, prof, alts)
N = calc_chapman_prof(alts, max(prof), alts(prof == max(prof)), H, 1);
err = sum((N - prof).^2);
% fprintf('%1.1f, %1.1e err \n', H, err)

