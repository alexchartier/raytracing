function D = load_nc(fn)

%%
% fn = '/Users/chartat1/data/superdarn/grid/20151220.bks.v3.0.grid.nc';

ni = ncinfo(fn);
D = [];
for i = 1:length(ni.Variables)

    varname = ni.Variables(i).Name;
    fname = strrep(varname, '.', '_');
    
    D(1).(fname) = double(ncread(fn, ni.Variables(i).Name));
end