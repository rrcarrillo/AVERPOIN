% [cone_k,cone_v]=trian_region(cone_elem) calculates the facets
% corresponding to the intersection polytope (except those that contain
% the origin vertex) by accessing the intersection points stored in
% cone_elem. These facets can be used as a decomposition of the region into
% simplices.
% cone_k and cone_v will contain the resulting facets and v-representation
% respectively.
function [cone_k,cone_v]=trian_region(cone_elem)
cone_ndims=length(cone_elem)-1; % Dimensions of the space that the cone spans
inter_v=[];
for ndim=1:cone_ndims % For each type of cone element: vertices, edges, faces,...
   nelems=length(cone_elem{ndim+1}); % Num. of cone elements of this dimensionality
   for nelem=1:nelems % For each element of this type
      inter_v=[inter_v cone_elem{ndim+1}(nelem).intersections]; % Collect all the intersection points
   end
end
% truncate and remove duplicated intersection vertices
inter_v(inter_v>1)=1;
inter_v(inter_v<0)=0;
cone_v=unique_tol(inter_v);

% Calculate the convex hull of the integration region
% Default option are {'Qt','Qx'} but these options (although generate a
% reduced number of facets) do not always work.
% Option {'QJ'} should be used to ensure a clearly-convex output although
% it is slower.
% 'QJ' option joggles each input coordinate by adding a random number to avoid precision errors.
% 'Qt' default option triangulates all non-simplicial facets before generating results.
% 'Qx' default option (n>4) merges a point into a coplanar facet, merges concave facets, merges duplicate ridges, and merges flipped facets.
% 'Pp' option does not report precision problems (warnings). It removes the narrow hull warning.
cone_k_hdim=convhulln0(cone_v,{'QJ','Pp'});

% num. of dimensions of the space that the cone spans, that is, num. of
% points per facets
red_cone_ndims=size(cone_k_hdim,1);

% Find all the facets which do not contain the origin vertex
cone_k_hdim_non_common=cone_k_hdim(:,all(cone_k_hdim~=0,1));
cone_nfacets=size(cone_k_hdim_non_common,2);
if cone_nfacets==0 % No proper facets of region cone have been found: degenerated input polytope?
   cone_k_vertices=unique(cone_k_hdim); % Probably too many (coplanar) vertices will be included
   cone_k=cone_k_vertices(cone_k_vertices~=0); % Remove origin
else
   cone_k=reshape(cone_k_hdim_non_common,red_cone_ndims,cone_nfacets); % recover shape after remove
end
