function d = dist(txloc, rxloc)

SphV = [(txloc(3) + 6371) * 1E3, deg2rad(txloc(1)), deg2rad(txloc(2)); 
        (rxloc(3) + 6371) * 1E3, deg2rad(rxloc(1)), deg2rad(rxloc(2))];
CartV = sphcart(SphV);

d = sqrt(sum((CartV(1, :) - CartV(2, :)) .^ 2)) / 1E3;

