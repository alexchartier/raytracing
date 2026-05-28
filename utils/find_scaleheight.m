function H = find_scaleheight(prof, alts, topside)
%% H = find_scaleheight(prof, alts, topside)
% Solves for the topside or bottomside scale height
% Just give it the portion of the profile you care about

%% optimize
% minimize (truth-model) dist to rx
% NOTE: I tried adding a stop condition, but it slows down the optimization
options = optimset('MaxFunEvals', 1000, 'TolX', 1, 'Display', 'On');

x0 = 100;
hmf2 = max(alts(prof == max(prof)));
if topside
    alti = alts >= hmf2;
else
    alti = alts <= hmf2;
end
f = @(X)iono_err(X, prof(alti), alts(alti));
[H, fval] = fminsearch(f, x0); % , options);

% disp(1)
% 
% N = calc_chapman_prof(alts, max(topside_prof), hmf2, X, 1);
% hold on
% plot(N, alts, 'r');
% plot(prof, alts, 'k');
% hold off

end

