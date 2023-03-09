%% load file
fname = 'data/wspr/2019-Mar-15-2200.nc';
data = read_netcdf(fname)

%% Plotting
rotang = -30;

txloc = links{2}.txlocs;
rxloc = links{2}.rxlocs;
home = links{2}.home;

latI_bad = [];
lonI_bad = [];
latI_good = [];
lonI_good = [];
for l = 1:length(home)
    % Great circle calcs
    [latI,lonI] = interpm([txloc(l, 1), rxloc(l, 1)], ...
        [txloc(l, 2), rxloc(l, 2)], 1, 'gc');
    if home(l) == 0
        latI_bad = [latI_bad; latI; NaN];
        lonI_bad = [lonI_bad; lonI; NaN];
    else
        latI_good = [latI_good; latI; NaN];
        lonI_good = [lonI_good; lonI; NaN];
    end
end

rotang = -30;
earth_example
% Plot the links
CartL = sphcart([64E5 * ones(size(latI_bad)), deg2rad(latI_bad), deg2rad(lonI_bad)]);
h3 = plot3(CartL(:, 1), CartL(:, 2), CartL(:, 3), '-r', 'LineWidth', 3);
rotate(h3, [0 0 1], rotang);


CartL = sphcart([64E5 * ones(size(latI_good)), deg2rad(latI_good), deg2rad(lonI_good)]);
h3 = plot3(CartL(:, 1), CartL(:, 2), CartL(:, 3), '-g', 'LineWidth', 3);
rotate(h3, [0 0 1], rotang);
%CartL = sphcart([64E5 * ones(size(latI_f)), deg2rad(latI_f), deg2rad(lonI_f)]);

rotate(globe, [0 0 1], rotang)