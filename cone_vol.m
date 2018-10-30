% vol=cone_vol(cone_elem) calculates the volume of the intersection polytope of cone_elem.
% The cone polytope is encoded by cell arrays with the format:
% elem{ndim+1}(nface).vertices which contains the vertices indices of the
% element nface of dimension ndim, except elem{1}(nvertex) which contains
% the vertex coordinate matrix.
% The original cone as explained in devel_cone().
function vol=cone_vol(cone_elem)
cone_ndims=length(cone_elem)-1; % Dimensions of the space that the cone spans
if cone_ndims > 1 && length(cone_elem{2}) > 0 % The cone must have at least one vertex
   space_ndims=size(cone_elem{1},1); % number of coordinates of the cone vertices (number of space dimensions in which the cone is)
   
   [cone_k,cone_v]=trian_region(cone_elem); % Decompose into simplices
   
   [n_vertices_facet,n_simplices]=size(cone_k); % number of simplicial cones in which the cone is divided
   vol=0;
   if n_vertices_facet == space_ndims % full-dimensional region polytope found
      for n_simplex=1:n_simplices % for each simplicial cone
         simplex_v=cone_v(:,cone_k(:,n_simplex));
         if rank(simplex_v)<space_ndims
            simplex_oriented_volume=0; % non-full-dimensional simplex
         else
            simplex_oriented_volume=det(simplex_v);
         end
         vol=vol+abs(simplex_oriented_volume);
      end
      vol=vol/factorial(space_ndims); % multiply by M!/(M+d)!
   end % if the polytope is non full-dimensional, the volume is 0
else
   if cone_ndims > 0 % In 1D
      vol=~isempty(cone_elem{2}); % Volume=1 if the cone has any ray (different from 0)
   else
      vol=0; % In 0D assume 0 volume
   end
end
