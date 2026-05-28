
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

%% 
n_O = sqrt( 1 - X ./ ...
    (1 - 1/2 .* Y.^2 .* sin(theta).^2 ./ (1 - X) + ...
    1 ./ (1 - X) .* sqrt(1 / 4 .* Y.^4 .* sin(theta).^4 + Y.^2 .* cos(theta).^2 .* (1 - X).^2) ...
    ));

n_X = sqrt( 1 - X ./ ...
    (1 - 1/2 .* Y.^2 .* sin(theta).^2 ./ (1 - X) - ...
    1 ./ (1 - X) .* sqrt(1 / 4 * Y.^4 .* sin(theta).^4 + Y.^2 .* cos(theta).^2 .* (1 - X).^2) ...
    ));


%% at the reflection point, N = 0

syms omega omega_0 omega_H theta
solve(1 ==  (omega_0 .^ 2) ./ (omega .^ 2) ./ ...
    (1 - 1/2 .* omega_H ./ omega.^2 .* sin(theta).^2 ./ (1 - (omega_0 .^ 2) ./ (omega .^ 2)) + ...
    1 ./ (1 - (omega_0 .^ 2) ./ (omega .^ 2)) .* sqrt(1 / 4 .* omega_H ./ omega.^4 .* sin(theta).^4 + omega_H ./ omega.^2 .* cos(theta).^2 .* (1 - (omega_0 .^ 2) ./ (omega .^ 2)).^2) ...
    ), omega)
