function prof = interp_sami(sami, lati, loni)
%% 3D stepwise interpolation of SAMI3 density to get a vertical profile

%TODO make sami grid bigger here or elsewhere
if loni < 0
    loni = loni + 360;
end

prof = zeros(size(sami.alt)) * NaN;
for ai = 1:length(sami.alt)
    prof(ai) = interp2(sami.lat, sami.lon, squeeze(sami.dene(ai, :, :))', lati, loni);
end
