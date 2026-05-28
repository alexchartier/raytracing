function data = read_netcdf(fname)

%% read generic netCDF into a struct
info = ncinfo(fname);
vn = {info.Variables(:).Name};
data = [];
for i = 1:length(vn)
    varname = vn{i};
    data.(varname) = ncread(fname, varname);
end



