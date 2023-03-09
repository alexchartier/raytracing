%% pharlap stuff
pl_path = '/Users/chartat1/pharlap/';
path(pl_path,path);
path([pl_path, 'src/matlab/'],path);
path([pl_path, 'mex/'],path);
setenv('DIR_MODELS_REF_DAT', [pl_path, 'dat/'])

% Set path print
fprintf('_______________________________________________________________________\n\n');
fprintf(' Set path to PHaRLAP: %s\n', pl_path);
fprintf('_______________________________________________________________________\n\n');



%% m_map
path(['/Users/chartat1/midas-3/matlab/m_map'],path);