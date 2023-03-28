function rms_err = rmse(errs)

bias = mean(errs);
rms_err = sqrt(mean((errs - bias) .^2));