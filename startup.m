%% pharlap stuff

[~, username] = system('whoami');
if isunix & ~ismac % must be on the Linux cluster
    pl_path = '/disks/d0378/users/cantrce1/Work/Models/pharlap_4.7.0/';
    base = sprintf('/disks/d0378/users/%s/data/nebula/', username);
else

    if strcmp(strip(username), 'chartat1')
        pl_path = '/Users/chartat1/pharlap/';
        path('~/MBeautifier/', path);
        base = '/Users/chartat1/data/nebula/';
    else
        pl_path = '/Users/cantrce1/Desktop/Work/Models/pharlap_4.7.0/';
        base = '/Users/cantrce1/Desktop/Work/Projects/Nebula/';
    end
end
path(pl_path, path);
addpath([pl_path, 'mex/'], '-begin')
addpath([pl_path, 'src/matlab/'], '-end')
setenv('DIR_MODELS_REF_DAT', [pl_path, 'dat/'])

path('./utils/', path);
path('./examples/', path);
path('./graphics/', path)

% Set path print
fprintf('_______________________________________________________________________\n\n');
fprintf(' Set path to PHaRLAP: %s\n', pl_path);
fprintf('_______________________________________________________________________\n\n');


% Plot
set(groot, ...
    'DefaultAxesFontSize', 20, ...
    'DefaultTextFontSize', 20, ...
    'DefaultTextFontName', 'Futura', ...
    'DefaultAxesFontName', 'Futura', ...
    'defaultfigurecolor', [1, 1, 1])
