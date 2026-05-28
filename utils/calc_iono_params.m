function [nmf2, hmf2, H] = calc_iono_params(n, alt)
%% calc_iono_params
% Calculate peak density, peak height and scale height from Ne profile
% % [nmf2, hmf2, H] = calc_iono_params(n, alt)
% plot(n, alt)
% plot(nmf2, hmf2, 'rx')
% text(nmf2, 500, sprintf('H: %1.0f km', H))

%% get params

[nmf2, I] = max(n);
hmf2 = alt(I);
H = interp1(n(I:end), alt(I:end) - hmf2, nmf2(1) / exp(1));
