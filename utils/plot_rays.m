function plot_rays(rays, varargin) %, txloc, rxloc) %, rot, ang)
%%
% plot_rays(rays, txloc, rxloc); % or
% plot_rays(rays{i}{:}, rays{i}{1}.txloc, rays{i}{1}.rxloc);
% earth_example
% rotate(gl, rot, ang)

%% Set optimization options
color = 0;
if length(varargin) == 1
    txloc = varargin{1};
    rxloc = txloc;
        color = 'r';

elseif length(varargin) == 2
    txloc = varargin{1};
    rxloc = varargin{2};
        color = 'r';

elseif length(varargin) == 3
    txloc = varargin{1};
    rxloc = varargin{2};
    color = varargin{3};
elseif length(varargin) == 0
    txloc = [0, 0, 0];
    rxloc = [0, 0, 0];
    color = 'r';
end

if iscell(rays)
    for r = 1:length(rays)
        if length(varargin) > 0
            plot_rays(rays{r}, txloc, rxloc, color);
        else
            plot_rays(rays{r});
        end
            
    end
    return
end

hold on
Re = 6371E3;
for r = 1:length(rays)
    if ~isempty(rays(r).initial_elev)
        hidx = rays(r).height >= -10;
        sph = [rays(r).height(hidx) * 1E3 + Re; deg2rad(rays(r).lat(hidx)); deg2rad(rays(r).lon(hidx))];
        cart = sphcart(sph');
        h3 = plot3(cart(:, 1), cart(:, 2), cart(:, 3), color);

    end
end

if length(varargin) > 0

    cart = sphcart([txloc(3) * 1E3 + Re, deg2rad(txloc(1)), deg2rad(txloc(2))]);
    h3 = plot3(cart(1), cart(2), cart(3), 'ro', 'markersize', 10, 'markerfacecolor', 'g');

    cart = sphcart([rxloc(3) * 1E3 + Re, deg2rad(rxloc(1)), deg2rad(rxloc(2))]);
    h3 = plot3(cart(1), cart(2), cart(3), 'ro', 'markersize', 10, 'markerfacecolor', 'r');

end
hold off





