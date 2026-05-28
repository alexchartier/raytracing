function n = appleton_hartree(OX_mode, Ne, B, wavevec, wavefreq)

% OX_mode =  0 for o, 1 for X
% theta = angle between wavevec and B

m_e = 9.10938356E-31; 
charge_e = 1.60217662E-19;
e_0 = 8.85418782E-12; 
plasmafreq = sqrt( (Ne * charge_e .^ 2) / (e_0 .* m_e.^2) ); 
omega = 2 * pi * wavefreq;  % angular freq of wave
omega_0 = 2 * pi * plasmafreq;  % angular plasma freq
omega_H = abs(B) * abs(charge_e) / m_e;   % electron gyrofreq


X = (omega_0 .^ 2) / (omega .^ 2); 
Y = omega_H / omega; 


%% Note not sure how to determine which is which
n_O = sqrt( (1 - X * (1 - X)) / ...
    (1 - X - 1/2 * Y^2 * sin(theta)^2 + ...
    sqrt((1 / 2 * Y^2 * sin(theta)^2)^2  + (1 - X)^2 * Y^2 * cos(theta)^2) ...
    ) ...
    );

n_X = sqrt( (1 - X * (1 - X)) / ...
    (1 - X - 1/2 * Y^2 * sin(theta)^2 - ...
    sqrt((1 / 2 * Y^2 * sin(theta)^2)^2  + (1 - X)^2 * Y^2 * cos(theta)^2) ...
    ) ...
    );