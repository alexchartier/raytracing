function [freqs, rg] = calc_ionogram(homed_rays)
% plot_ionogram(homed_rays)

arguments 
    homed_rays cell
end

nfreq = length(homed_rays);
maxnrg = 0;
for ri = 1:length(homed_rays)
    rays = homed_rays{ri};
    if length(rays) > maxnrg
        maxnrg = length(rays);
    end
end

freqs = zeros(nfreq, 1);
rg = zeros(nfreq, maxnrg) * NaN;
for ri = 1:length(homed_rays)
    rays = homed_rays{ri};
    freqs(ri) = rays(1).frequency;
    
    for r = 1:length(rays)
        rg(ri, r) = rays(r).group_range_to_rx;
    end
end