%% Remap SAMI3 to a regular geographic grid - this is for the HIAMCM version
% note alt is the second dimension
%in_fname = 'data/sami/sami3_mar_2015.nc';
clear

%% Filenames
in_fn_fmt = '/Users/chartat1/data/sami3/2017_tid/sami3_{ddmmmYYYY}.nc';
out_fn_fmt = '/Users/chartat1/data/sami3/2017_tid/sami_mat/{YYYY-mm-dd_HHMM}.mat';
stimes = datenum(2017, 1, 10:12);

%% Loop over days
for d = 1:length(stimes)
    stime = stimes(d);

    %% load
    D = load_sami(filename(in_fn_fmt, stime), {'time', 'dene0', 'lat0', 'lon0', 'alt0'}, stime);

    %% Setup storage
    D2.lat = -88:2:88;
    D2.lon = 0:5:360;
    D2.alt = 92:2:800;
    D.time = round(D.time * 60*24) / (60*24);

    %% Loop and interpolate
    for t = 1:length(D.time)
        %% Zeros
        D2.time = D.time(t);
        D2.dene = zeros(length(D2.alt), length(D2.lat), length(D2.lon));

        %% Loop over alts to get a regular lat/lon
        dene_t = zeros(length(D.alt0), length(D2.lat), length(D2.lon));
        for a = 1:length(D.alt0)
            lat = D.lat0(:, a, :);
            lon = D.lon0(:, a, :);
            dene_in = D.dene0(:, a, :, t);
            dene_t(a, :, :) = griddata(lat(:), lon(:), dene_in(:), D2.lat, D2.lon')';
        end

        % clean up dateline
        dateline = (dene_t(:, :, 2) + dene_t(:, :, end-1)) ./ 2;
        dene_t(:, :, 1) = dateline; 
        dene_t(:, :, end) = dateline;
        
        %% Loop over lat/lon to get a regular alt
        for l1 = 1:length(D2.lat)
            for l2 = 1:length(D2.lon)
                D2.dene(:, l1, l2) = interp1(D.alt0, dene_t(:, l1, l2), D2.alt);
            end
        end

        %% Calculate PHARLAP stuff
        [D2.iono_en_grid, ~, D2.collision_freq, D2.iono_grid_parms, ...
            D2.Bx, D2.By, D2.Bz, D2.geomag_grid_parms] = gen_grid_parms(D2);
        
        %% Save
        out_fn = filename(out_fn_fmt, D2.time);
        savestruct(out_fn, D2)
        fprintf('Saved to %s\n', out_fn)

    end
    clear D D2 dene_t dene_in lat lon

end
%
%
% %% test scatteredInterpolant
% F = scatteredInterpolant(D.lat(:), D.lon(:), nmf2_t(:), 'natural', 'linear');
% [lon, lat] = meshgrid(D2.lon, D2.lat);
% nmf2 = F(lat(:), lon(:));
% nmf2 = reshape(nmf2, length(D2.lat), length(D2.lon));
% contourf(nmf2)
%


















