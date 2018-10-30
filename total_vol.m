% vol=total_vol(cone_elem) calculates the volume of each intersection polytope
% stored in cone_elem.
% vol is a vector containing the volume of the cone and all the regions
% cones intersected with the hypercube.
% The cone polytope is encoded by cell arrays with the format:
% elem{ndim+1}(nface).vertices which contains the vertices indices of the
% element nface of dimension ndim, except elem{1} which contains
% the vertex coordinate matrix.
function vol=total_vol(cone_elem)
cone_ndims=length(cone_elem)-1; % Dimensions of the space that the cone spans
if cone_ndims==0 % All the rays are zero
   vol=1; % Volume of the hypercube
else
   cone_vol_val=cone_vol(cone_elem);
   vol=[cone_vol_val]; % The volume of each region will be appended in this array
   for ndim=1:(cone_ndims-1) % For each type of cone element: vertices, edges, faces,...
      nelems=length(cone_elem{ndim+1}); % Num. of cone elements of this dimensionality
      for nelem=1:nelems % For each element of this type
         reg_vol=cone_vol(cone_elem{ndim+1}(nelem).region);
         vol=[vol reg_vol];
      end
   end
end
