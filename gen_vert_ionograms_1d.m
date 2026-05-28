in_fn_fmt = '/Users/chartat1/data/sami3/2017_tid/sami_mat/{YYYY-mm-dd_HHMM}.mat';
out_fn_fmt = '/Users/chartat1/data/sami3/gs_ionograms/ionogram_%i.mat';
time = datenum(2017, 1, 10, 18, 0, 0);
OX_mode = 1;
freqs = 2:0.1:15;
txalt = 580;

%% Load
sami = loadstruct(filename(in_fn_fmt, time));
alts = sami.alt;
%% Loop through and calculate ionograms
ct = 0;
for l1 = 1:length(sami.lat)
    lat = sami.lat(l1);
    for l2 = 1:length(sami.lon)
        lon = sami.lon(l2);

        %% mag field and Ne profile
        [~, ~, ~, Inc, B] = igrfmagm(alts, lat * ones(size(alts)), ...
            lon * ones(size(alts)), str2num(datestr(time, 'YYYY')) * ones(size(alts))); 
        Inc = Inc';  % inclination
        B = B' / 1E9;  % field strength

        Ne = sami.dene(:, l1, l2) * 1E6;
        %% Skip if electron density profile is bad
        if sum(isnan(Ne)) > 0
            continue
        end

        %% raytrace
        vht = raytrace_1d(alts, freqs, Ne, B, Inc, OX_mode, txalt);


        %% Store
        % out.ionogram = sparse(length(freqs), length(vht));
        % for f = 1:length(freqs)
        %      out.ionogram(f, vht == round(rg(f))) = 1;
        % end
        out.freqs = freqs;
        out.vht = vht; 
        out.profile = Ne;
        out.alts = alts;
        out.txalt = txalt;
        out.nmf2 = max(Ne);
        out.hmf2 = alts(Ne == max(Ne));
        savestruct(sprintf(out_fn_fmt, ct), out)
        fprintf('Saved to %s\n', sprintf(out_fn_fmt, ct))
        ct = ct + 1;
    end
end
