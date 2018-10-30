% cone_elem=region_intersec(cone_elem,cube_elem) calculates the intersections
% of all the cone elements (cone, facets, ridges, peaks,..., faces, edges and
% vertices) for the cones corresponding to the regions of cone_elem with 
% their correponding hypercube element (vertices, edges,..., facets).
% All the polytopes are encoded by cell arrays with the format:
% elem{ndim+1}(nface).vertices which contains the vertices indices of the
% element nface of dimension ndim, except elem{1} which contains the
% vertex coordinate matrix.
% The cone is encoded by the polytope resulting from its projection on the
% hyperplane with director vector inside the cone.
% Intersection point are stored in the "intersection" field of each
% element.
% The region cones are expected to be found in the "region" cone of each
% input cone element
function cone_elem=region_intersec(cone_elem,cube_elem)
cone_ndims=length(cone_elem)-1; % Num. of dimensions of the space that the cone spans
for cone_ndim=1:(cone_ndims-1) % For each type of cone element: rays, faces,..., except cone (cone_ndims-1)
   % In the cone, cone_elem{2} are vertices (rays) (dim. 1)
   cone_nelems=length(cone_elem{cone_ndim+1}); % Num. of cone elems. of this dim.
   for cone_nelem=1:cone_nelems % For each element of this type
      % Calculate the intersection of the current region cone with the hypercube
      cone_elem{cone_ndim+1}(cone_nelem).region=cone_intersec(cone_elem{cone_ndim+1}(cone_nelem).region,cube_elem);
   end   
end
