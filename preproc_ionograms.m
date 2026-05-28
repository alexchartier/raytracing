%% preproc_ionograms.m
% Convert the input dataset to a smaller format, and clean up
% Keep the el and az
% See also gen_vert_ionograms.m, gen_oblique_ionograms.m, correlate_ionograms.m

%% Set inputs
times = datenum(2017, 1, 10, 12, 0, 0):3/24:datenum(2017, 1, 13, 0, 0, 0);

vert = true;
if vert
    in_dir_fmt = '/Users/chartat1/data/sami3/gs_ionograms_3d/hiamcm/{YYYY-mm-dd_HHMM}/vert/';
    out_fn = '/Users/chartat1/data/sami3/gs_ionogram_set_vert.mat';
else
    in_dir_fmt = '/Users/chartat1/data/sami3/gs_ionograms_3d/hiamcm/{YYYY-mm-dd_HHMM}/oblique/';
    out_fn = '/Users/chartat1/data/sami3/gs_ionogram_set_oblique_600.mat';
end
mod_fn_fmt = '/Users/chartat1/data/sami3/2017_tid/sami_mat/{YYYY-mm-dd_HHMM}.mat';

sl_cutoff_pct = 10;


%% Loop over times
iono_set = {};
ct = 1;

BS.H_bs = [];
BS.hmax = [];
BS.nmax = [];
BS.lat = [];
BS.lon = [];
BS.lt = [];


for t = 1:length(times)

    %% Generate ionogram file-list and load model for that time
    flist = dir(filename(in_dir_fmt, times(t)));
    flist = flist(3:end);

    sami = loadstruct(filename(mod_fn_fmt, times(t)));

    %% Loop over files and store good ionograms
    for fi = 1:length(flist)
        fname = [flist(fi).folder, '/', flist(fi).name];
        arr = strsplit(flist(fi).folder, '/');
        rays = loadstruct(fname);
        iono = get_ionogram(rays, sl_cutoff_pct);
        iono.time = datenum(datetime(arr{8},'InputFormat', 'yyyy-MM-dd_HHmm'));

        if length(iono.freqs(:)) < 10 % skip empty
            disp('empty')
            continue
        end

        iono.mod_prof = interp_sami(sami, iono.lat_ref, iono.lon_ref);
        iono.prof_alts = sami.alt;
        [iono.nmax, iono.fmax, iono.hmin, iono.grouprg, ray] = id_golden_ray(rays);
        iono.H_bs = find_scaleheight(iono.mod_prof, sami.alt, false);
        [~, maxi] = max(iono.mod_prof);
        iono.hmax = sami.alt(maxi);
        iono_set{ct} = iono;

        ct = ct + 1;

    end

    %% Grab the bottomside scale height params across the whole model
    for i = 1:length(sami.lat)
        for j = 1:length(sami.lon)
            prof = sami.dene(:, i, j);
            if sum(isnan(prof)) > 0
                continue
            end
            H_bs = find_scaleheight(prof', sami.alt, false);
            [nmax, maxi] = max(prof);
            hmax = sami.alt(maxi);
            lat = sami.lat(i);
            lon = sami.lon(j);
            lt = ((time - floor(time)) + lon / 360) * 24;
            lt(lt < 0) = lt(lt < 0) + 24;
            lt(lt >= 24) = lt(lt >= 24) - 24;

            BS.H_bs = [BS.H_bs, H_bs ];
            BS.hmax = [BS.hmax, hmax];
            BS.nmax = [BS.nmax, nmax];
            BS.lat = [BS.lat, lat];
            BS.lon = [BS.lon, lon];
            BS.lt = [BS.lt, lt];

        end
    end

end

%% Calculate relations between fmax and nmax, grouprg and hmin, hmax and H_bs
fmax = zeros(size(iono_set));
nmax = zeros(size(iono_set));
hmin = zeros(size(iono_set)); % the lowest ray height
hmax = zeros(size(iono_set));
H_bs = zeros(size(iono_set));
grouprg = zeros(size(iono_set));
lat = zeros(size(iono_set));
lt = zeros(size(iono_set));


for i = 1:length(iono_set)
    fmax(i) = iono_set{i}.fmax;
    nmax(i) = iono_set{i}.nmax;
    hmin(i) = iono_set{i}.hmin;
    hmax(i) = iono_set{i}.hmin;
    H_bs(i) = iono_set{i}.H_bs;
    grouprg(i) = iono_set{i}.grouprg;
        lat(i) = iono_set{i}.lat_ref;
              lt(i) = ((time - floor(time)) + iono_set{i}.lon_ref / 360) * 24;
end

lt(lt < 0) = lt(lt < 0) + 24;
lt(lt >= 24) = lt(lt >= 24) - 24;
iono_set{1}.mdl_fmax_sqrtnmax = fitlm(fmax, sqrt(nmax));
iono_set{1}.mdl_grouprg_peakht = fitlm(grouprg, hmin);
iono_set{1}.mdl_hmax_H_bs = fitlm(hmax, H_bs);    

% X =  [ones(size(BS.H_bs)); BS.lat; BS.lon; BS.hmax; BS.nmax; BS.lt];
ltmod = BS.lt - 15;
ltmod(ltmod < 0) = ltmod(ltmod < 0) + 15;
X =  [abs(BS.lat); BS.hmax; ltmod];
iono_set{1}.mdl_Hbs_abslat_hmax_ltmod = regress(BS.H_bs', X');

%%  Save
savestruct(out_fn, iono_set)
fprintf('Saved to %s\n', out_fn)




















