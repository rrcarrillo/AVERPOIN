% cone_elem=add_regions(cone_elem) for each cone element (vertices, edges, ...
% facets) calculates a cone which defines the space normal to that element
% and that is not shared by other elements. That is, the space of points that
% when proyected on the element are only inside that element.
% The d rays of the cone defining this space are the original cone k element
% rays plus the (d-k) vectors normal to the adjacent element facets.
% The generated cone corresponding to each region is encoded into a data
% structure simular to cone_elem and inserted in each original cone
% element.
% The cone polytope is encoded by cell arrays with the format:
% elem{ndim+1}(nface).vertices which contains the vertices indices of the
% element nface of dimension ndim, except elem{1}(nvertex) which contains
% the vertex coordinate matrix.
% The original cone is encoded by the polytope resulting from its projection
% on a hyperplane with director vector the cone center ray.
function cone_elem=add_regions(cone_elem)
cone_ndims=length(cone_elem)-2; % Dimensions of the space that the (projected) cone spans

if cone_ndims > 0 && ~isempty(cone_elem{2}) % The cone must have at least one vertex
   space_ndims=size(cone_elem{1},1); % Check the number of coordinates of the cone vertices to find out the number of space dimensions in which the cone is
   % Since the cone is encoded as a projection of the cone on a hyperplane,
   % the number of space dimensions must be equal to the number of cone
   % element types plus 1
   if space_ndims ~= cone_ndims+1
      error('Wrong cone: Unexpected number of space dimension according to the number of element types')
   end % Consider degenerate cases?
   
   % Include normal vector in each facet element
   cone_elem=add_normals(cone_elem);

   for ndim=1:cone_ndims % For each type of cone element: vertices, edges, faces,...
      % In the cone, cone_elem{1} is the vertex coordinate matrix referenced by all the cone elements
      nelems=length(cone_elem{ndim+1}); % Num. of cone elements of this dimensionality
      for nelem=1:nelems % For each element of this type
         ind_elem_vertices=cone_elem{ndim+1}(nelem).vertices;
         ind_elem_facets=cone_elem{ndim+1}(nelem).facets;
         coord_elem_vertices=cone_elem{1}(:,ind_elem_vertices); % vertex coordinates (rays) of the current element
         coord_elem_normals=[cone_elem{end-1}(ind_elem_facets).normal]; % coordinates of the vectors which are normal to the facets which contain the current element
         elem_cone_v=[coord_elem_vertices coord_elem_normals]; % rays of the region cone
         
         % Create the data structure for the cone defining the region and
         % store it in the "region" field of the corresponding cone element 
         cone_elem{ndim+1}(nelem).region=devel_cone(elem_cone_v);
      end
   end
end

