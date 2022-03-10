function Cart = sphcart( Sph, SphV )
% Function to convert spherical to cartesian coordinates
%    [Cart,CartV] = sphcart( Sph, SphV )
% Input Arguments:
%    Sph  = Spherical coordinates, [rad(m),lat(rad),lon(rad);...]
% Optional Input Arguments:
%    SphV = Spherical vector, [Up(m),North(rad),East(rad)]
% Output arguments:
%    Cart = Cartesian coordinates, [x,y,z;...]
% Optional Output Arguments:
%    CartV = Cartesian vector, [vx,vy,vz;...]   
% Notes:
%    For geographic coordinates, x points to Greenwich, z to the north pole
% Coordinates or vectors can be passed as cell arrays {x,y,z}. If Sph is a
% cell array the result is returned in a cell array.
%
% See also CARTSPH
% -------------------------------------------------------------------------
if iscell(Sph)
   Cart = { Sph{1}.*cos(Sph{2}).*cos(Sph{3}),...
            Sph{1}.*cos(Sph{2}).*sin(Sph{3}),...
            Sph{1}.*sin(Sph{2}) };
else
   Cart = [ Sph(:,1).*cos(Sph(:,2)).*cos(Sph(:,3)),...
            Sph(:,1).*cos(Sph(:,2)).*sin(Sph(:,3)),...
            Sph(:,1).*sin(Sph(:,2)) ];
end
if nargin == 1, CartV = []; return; end

% Vector transformation
if iscell(Sph)
   if ~iscell(SphV) 
      SphV = {reshape(SphV(:,1),size(Sph{1})),reshape(SphV(:,2),size(Sph{1})),...
              reshape(SphV(:,3),size(Sph{1}))};
   end
   Cart = repmat({SphV{1}},1,3).*...
             {cos(Sph{3}).*cos(Sph{2}),sin(Sph{3}).*cos(Sph{2}),sin(Sph{2})}+...
          repmat({SphV{2}},1,3).*...
             {-cos(Sph{3}).*sin(Sph{2}),-sin(Sph{3}).*sin(Sph{2}),cos(Sph{2})}+...
          repmat({SphV{3}},1,3).*...
             {-sin(Sph{3}),cos(Sph{3}),zeros(size(Sph{1}))};
else
   Cart = repmat(SphV(:,1),1,3).*...
             [cos(Sph(:,3)).*cos(Sph(:,2)),sin(Sph(:,3)).*cos(Sph(:,2)),sin(Sph(:,2))]+...
          repmat(SphV(:,2),1,3).*...
             [-cos(Sph(:,3)).*sin(Sph(:,2)),-sin(Sph(:,3)).*sin(Sph(:,2)),cos(Sph(:,2))]+...
          repmat(SphV(:,3),1,3).*...
             [-sin(Sph(:,3)),cos(Sph(:,3)),zeros(length(Sph(:,3)),1)];
end
