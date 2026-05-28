function params = fit_curve(x, y)
%% fit_curve.m
% Fit a hyperbola to a set of points
% params = fit_curve(x, y)
% plot(x, y, 'pg')
% hold on
% plot(x, hyperbola(params, x), '-r')
% hold off
% grid

xoff = min(x);
yoff = min(y);
B0 = [1, 1, 1];

NRCF = @(b) norm((y - yoff)- hyperbola(b, x - xoff));   % Residual Norm Cost Function
mdl = fitnlm(x - xoff, y - yoff, @hyperbola, B0);  % fit nonlinear model

B = fminsearch(NRCF, B0);                               % Estimate Parameters

% plot(x - xoff, y - yoff, 'pg')
% hold on
% plot(x - xoff, hyperbola(B, x - xoff), '-r')
% hold off


params = [B(1), B(2) - xoff, B(3) + yoff];





