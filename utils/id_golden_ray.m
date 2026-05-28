function [nmax, fmax, hmin, grouprg, ray] = id_golden_ray(rays)

%% id_golden_ray
% Find ray that penetrates the deepest, get corresponding density

hmin = 1E9;
for r = 1:length(rays)
    for ri = 1:length(rays{r})
        ray = rays{r}(ri);
        if min(rays{r}(ri).height) < hmin
            hmin = min(rays{r}(ri).height);
            grouprg = rays{r}(ri).group_range_to_rx;
            id = rays{r}(ri).height == hmin;
            fmax = rays{r}(ri).frequency;
            nmax = rays{r}(ri).electron_density(id);
            ray = rays{r}(ri);
        end
    end
end