function N = calc_chapman_prof(h, Np, hp, H, c)
%% generate ionospheric profiles
% see Verhulst and Stankov (2013) for details
% they found exponential was the best for topside data, then chapman
% Exponential looks weird at peak though, so I'm going with Chapman
% 
% h - altitudes
% Np - peak density
% hp - peak height
% H - scale height
% c - chapman parameter (typically 0 - 1
% 
% Np = 1E12;
% hp = 301;
% H = 200;
% h = 150:2:600;
% c = 1; % chapman parameter
% N = calc_chapman_prof(h, Np, hp, H, c);
% plot(N, h, '--.')
% grid on
% grid minor
% 

% %% Set default args
% params = {'mode'};
% defaults = {'exp'};
% varargparse(varargin, params, defaults);

% %%  Exponential
% N = zeros(size(h));
% iht = h >= hp;
% Ntop = Np .* exp(- (h(iht) - hp) ./ H);
% Nbot = Np .* exp(- (hp - h(~iht)) ./ H);
% N(~iht) = Nbot;
% N(iht) = Ntop;

%% Chapman
% N = Np .* exp(c * (1 - (h - hp) ./ H - exp(- (h - hp) ./ H)));

prof = exp(c .* (1 - (h - hp) ./ H - exp(- (h - hp) ./ H)));

N = prof ./ max(prof) .* Np;


