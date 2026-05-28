%% subgrid sami function
function sami = subgrid_sami(sami, satlats, satlons)
sami.lon(sami.lon >= 180) = sami.lon(sami.lon >= 180) - 360;
lati = sami.lat >= min(satlats) - 5 & sami.lat <= max(satlats) + 5;
loni = sami.lon >= min(satlons) - 5 & sami.lon <= max(satlons) + 5;

sami.lon = sami.lon(loni);
sami.lat = sami.lat(lati);
sami.dene = sami.dene(:, lati, loni);

end