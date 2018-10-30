% cone_elem=add_normals(cone_elem) calculates a vector which is normal to
% each cone facets and points towards outside the cone. This vector is
% stored in the field "normal" in each cone facet.
% The cone polytope is encoded by a cell array with the format:
% elem{ndim+1}(nface).vertices which contains the vertices indices of the
% element nface of dimension ndim, except elem{1} which contains
% the coordinate matrix of the vertices.
% This cone is encoded by the polytope resulting from its projection on the
% hyperplane with director vector (1,1,...1) (See project_cone()).
function cone_elem=add_normals(cone_elem)
ndims=length(cone_elem)-2; % Num. of dimensions of the space which the cone spans

if ndims>0 % The polytope must have at least 1 dimensions
   nfacets=length(cone_elem{end-1}); % Number of cone facets
   
   % Calculate an approximate cone center ray
   cone_pcenter=mean(cone_elem{1},2);

   for nfacet=1:nfacets % For each cone facet
      % Find out the coordintates of the current-facet vertices (rays)
      facet_vertex_indices=cone_elem{end-1}(nfacet).vertices;
      facet_vertices=cone_elem{1}(:,facet_vertex_indices);

      normal_v=normal_vec(facet_vertices); % Vector normal (n) to the current facet

      outer_f=facet_orient(cone_pcenter,facet_vertices,normal_v); % Return true if normal_v points outside the cone

      % Store the normal vector in the facet element
      % The normal vector direction is inverted if it is initially pointing
      % towards inside the cone (outer_f==0)
      cone_elem{end-1}(nfacet).normal=normal_v*(-1+2*outer_f);
   end
end
