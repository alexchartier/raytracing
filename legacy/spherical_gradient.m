function [df_theta, df_phi] = spherical_gradient(rad, lat, lon, fn)
%% function to calculate the gradient in spherical coordinates (expect m, deg)


r = rad;
theta = deg2rad(90 - lat);
phi = deg2rad(lon);

[grad_theta, grad_phi] = gradient(fn, unique(theta), unique(phi));
df_theta = 1 ./ r .* grad_theta;
df_phi = 1 ./ (r .* sin(theta)) .* grad_phi;
