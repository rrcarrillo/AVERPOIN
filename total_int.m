% total_int=total_int(cone_elem) calculates the sum of the integrals of the
% squared distance function to the cone over each intersection polytope
% stored in cone_elem.
% The cone polytope is encoded by cell arrays with the format:
% elem{ndim+1}(nface).vertices which contains the vertices indices of the
% element nface of dimension ndim, except elem{1}(nvertex) which contains
% the vertex coordinate matrix.
function total_int=total_int(cone_elem)
cone_ndims=length(cone_elem)-1; % Dimensions of the space that the cone spans
if cone_ndims==0 % All the rays are zero
   total_int=size(cone_elem{1},1)/3; % Maximal integral value
else
   total_int=0; % Accumulated region integral values
   for ndim=1:(cone_ndims-1) % For each type of cone element: vertices, edges, faces,...
      nelems=length(cone_elem{ndim+1}); % Num. of cone elements of this dimensionality
      for nelem=1:nelems % For each element of this type
         subspace_v=cone_elem{1}(:,cone_elem{ndim+1}(nelem).vertices);
         reg_int=cone_int(cone_elem{ndim+1}(nelem).region,subspace_v);
         total_int=total_int+reg_int;
      end
   end
end
