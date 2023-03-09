function mix = load_mix(pot_fname)
%% Load_Mix.m 
% % function to load the MIX input
% pot_fname_fmt = '/Users/chartat1/pymix/data/pot_sami_cond/mar19/ampere_mix_{YYYY}-{mm}-{DD}T{HH}-{MM}-{SS}Z.nc';
% time = datenum(2019, 3, 2, 10, 0, 0);
% pot_fname = filename(pot_fname_fmt, time);


%%
mix = [];

keys = {'Potential', 'Geographic Latitude', 'Geographic Longitude', 'MLAT (AACGM)', 'MLON (AACGM)'};

for k = keys
    outname = regexprep(k{1}, ' ', '_');
    outname = regexprep(outname, '(', '');
    outname = regexprep(outname, ')', '');
    mix.(outname) = double(ncread(pot_fname, k{1}));
    mix.(outname) = mix.(outname)(2:end, :); 
end