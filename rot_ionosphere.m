%% rot_ionosphere.m
% Move the ionosphere in longitude

times = datenum(2017, 1, 5:6:23, 11, 0, 0);

verts = [0, 1];
indir = '/Users/chartat1/data/sami3/2017_tid/recon_iono/';

loninc = 5;
%%


flist = dir(indir);
flist = flist(3:end);


for f = 1:length(flist)
    in_fn = [indir, flist(f).name];
    out_fn = [indir, 'rot_', flist(f).name];
    D = loadstruct(in_fn); 
    D.lon = D.lon - loninc;
    D.iono_grid_parms(4) = D.iono_grid_parms(4) - loninc;
    D.geomag_grid_parms(4) = D.geomag_grid_parms(4) - loninc;
    savestruct(out_fn, D)
end