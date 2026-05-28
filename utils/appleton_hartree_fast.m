function n = appleton_hartree_fast(theta, Ne, B, wavefreq, OX_mode)
%% Calculate the complex refractive index from the Appleton-Hartree equation
%
% [n_O, n_X, n] = appleton_hartree(theta, Ne, B, wavefreq)
% 
% theta = angle between wavevec and B in radians
% Ne = electron density
% B = geomagnetic field strength
% wavefreq = EM wave frequency in Hertz

%% Constants
m_e = 9.10938356E-31; 
charge_e = 1.60217662E-19;
e_0 = 8.85418782E-12; 


%% Calculate parameters

omega = 2 * pi .* wavefreq;  % angular freq of wave
omega_0 = sqrt( (Ne .* charge_e .^ 2) / (e_0 .* m_e) );  % angular plasma freq
omega_H = abs(B) .* abs(charge_e) ./ m_e;   % electron gyrofreq

X = (omega_0 .^ 2) ./ (omega .^ 2); 
Y = omega_H ./ omega; 


%% Calculate refractive index

% O
if OX_mode == 1
n = sqrt( 1 - X ./ ...
    (1 - 1/2 .* Y.^2 .* sin(theta).^2 ./ (1 - X) + ...
    1 ./ (1 - X) .* sqrt(1 / 4 .* Y.^4 .* sin(theta).^4 + Y.^2 .* cos(theta).^2 .* (1 - X).^2) ...
    ));

% X
elseif OX_mode == -1
n = sqrt( 1 - X ./ ...
    (1 - 1/2 .* Y.^2 .* sin(theta).^2 ./ (1 - X) - ...
    1 ./ (1 - X) .* sqrt(1 / 4 * Y.^4 .* sin(theta).^4 + Y.^2 .* cos(theta).^2 .* (1 - X).^2) ...
    ));
else
% no B case
plasmafreq = omega_0 ./ (2 .* pi);
n = sqrt(1 - plasmafreq .^ 2 ./ wavefreq .^2);  
end
