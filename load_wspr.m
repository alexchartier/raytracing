function crd_o = load_wspr(fname)

%% wspr_txt comes from http://wsprnet.org/drupal/wsprnet/spotquery
% fname = 'wspr.txt';
% [lat_t, lon_t, lat_r, lon_r] = load_wspr(fname);
%%

txt = asciiread(fname);

lat_t = [];
lon_t = [];
lat_r = [];
lon_r = [];

for l = 1:size(txt, 1)
    crd = split(strip(txt(l, :)));
    if len(crd) > 1
        md_crd1 = upper(crd{1}(~isspace(crd{1})));
        md_crd2 = upper(crd{2}(~isspace(crd{2})));
        if length(md_crd1) == 6 && length(md_crd2) == 6
            [la_t, lo_t] = maidenhead_to_ll(md_crd1);
            [la_r, lo_r] = maidenhead_to_ll(md_crd2);
           
            if la_t > 40 && la_r > 40
                lat_t = [lat_t, la_t];
                lon_t = [lon_t, lo_t];
                lat_r = [lat_r, la_r];
                lon_r = [lon_r, lo_r];
            end
        end
    end 
end

%%
%crd_o = sort([lat_t; lon_t; lat_r; lon_r], 2, 'descend'); 
crd_o = [lat_t; lon_t; lat_r; lon_r]; 

