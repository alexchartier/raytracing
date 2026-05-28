function fit_rays = find_fit_rays(rays)

%% find_golden_rays
% Pick out the rays that are most important in characterizing the
% ionosphere
%
% TODO: deal with multi-trace cases


%% get the ionogram
freqs = [];
ranges = [];
for ri = 1:length(rays)
    if isstruct(rays{ri})
        r2 = rays{ri};
        for r = 1:length(r2)
            freqs = [freqs; r2(r).frequency];
            ranges = [ranges; r2(r).group_range_to_rx];
        end
    end
end

[ranges, I] = sort(ranges);
freqs = freqs(I);
new = diff(freqs);
freqs = freqs(new ~= 0);
ranges = ranges(new ~= 0);

%% #1 get the turning points
d1 = round(diff(freqs) * 10);
I = sign(d1(2:end)) ~= sign(d1(1:end-1));

ranges_s = ranges(2:end-1);
freqs_s = freqs(2:end-1);
% 
% hold on
% plot(ranges_s, freqs_s, '.k')
% plot(ranges_s(I), freqs_s(I), '.r', 'markersize', 30)
% hold off

%% Analyze and store
n_tpi = sum(I);  % determine number of turning points
flist = [];
for ri = 1:length(rays)
    flist = [flist; rays{ri}(1).frequency];
end

if all(rays{1}(1).txloc == rays{1}(1).rxloc)% -> VIS
    ris = flist == max(flist) | flist == max(flist) - 1;
    fit_rays = rays(ris);
    fit_rays = fit_rays(end:-1:1);  % flip to get the 'golden' ray in 1st position
elseif n_tpi == 0  % -> single-mode OIS
    ris = flist == min(flist) | flist == min(flist) + 2;
    fit_rays = rays(ris);
else % TODO -> classic oblique with low/high mode, or something else. need to figure out later
    min_rg = min(ranges_s(I));
end



























