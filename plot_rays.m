function plot_rays(txloc, rxloc, rays) %, rot, ang)
%%
% plot_rays(txloc, rxloc, rays)

% earth_example
% rotate(gl, rot, ang)

hold on
Re = 6371E3; 
for r = 1:length(rays)
    if ~isempty(rays(r).initial_elev)
    hidx = rays(r).height >= -10;
    sph = [rays(r).height(hidx) * 1E3 + Re; deg2rad(rays(r).lat(hidx)); deg2rad(rays(r).lon(hidx))];
    cart = sphcart(sph');
    h3 = plot3(cart(:, 1), cart(:, 2), cart(:, 3), 'm');
    %h4 = plot3(cart(:, 1), cart(:, 2), cart(:, 3), '.k');

    %rotate(h3, rot, ang)
    end
end

cart = sphcart([txloc(3) * 1E3 + Re, deg2rad(txloc(1)), deg2rad(txloc(2))]);
h3 = plot3(cart(1), cart(2), cart(3), 'ro', 'markersize', 10, 'markerfacecolor', 'r');
%rotate(h3, rot, ang);

cart = sphcart([rxloc(3) * 1E3 + Re, deg2rad(rxloc(1)), deg2rad(rxloc(2))]);
h3 = plot3(cart(1), cart(2), cart(3), 'go', 'markersize', 10, 'markerfacecolor', 'g');

%rotate(h3, rot, ang);

hold off





