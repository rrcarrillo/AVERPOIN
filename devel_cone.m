% elem=devel_cone(cone_rays)
% devel_cone generates all the elements (cones, facets, ridges, 
% peaks,..., faces and edges (rays) of the cone convex hull
% defined by cone_rays.
% cone_rays is a n x r matrix which cotains the coordinates of r rays in
% an n-dimensional space.
% devel_cone returns elem, which is an array of cells with the content:
% elem{1} is a n x r0 matrix which contains the ray coordinates of the
% cone convex hull.
% elem{ndim+1}(nelem).vertices is a ray-index column vector of length ndim
% which contains indices to the columns (rays) of elem{1} which define
% the element nelem of dimension ndim (ndim=1: rays, ndim=2: faces,...).
% (the origin vertex is not specified in the element definitions since it is
% assumed to be included in every cone element, thus elem will not contain
% the 0 vertex).
%
% To achieve this, this function performs these calculation:
% -The length of the input rays is truncated by the hyperplane whose
%  normal vector is near the cone center ray.
% -The convex hull of these points plus the coordinate origin is
% calculated, obtaining a convex polytope.
% -The facets of the cone convex hull from the polytope
% -The cone convex hull is divided into simplex cones
% -The rest of the cone convex hull elements (ridges, peaks, etc) are
%  calculated from the facets
function elem=devel_cone(cone_rays)
space_ndims=size(cone_rays,1); % Number of space dimensions
cone_rays=cone_rays(:,any(abs(cone_rays) > eps)); % Remove zero rays
cone_rays_uni=unique_tol(normalize_vecs(cone_rays)); % Remove duplicates
if size(cone_rays_uni,2) > 0 % We deal with the case in wich there are no input rays independently
   % Calculate the convex hull of the polytope defined by the cone rays plus the coord. origin vertex
   % The option {'QV0'} (plus options {'Qt','Qx'}) could be used to specify that only facets which
   % contain the origin vertex are good (returned), we save a few milliseconds and some code below
   conv_cone_k=convhulln0(cone_rays_uni,{'Qt','Qx','Pp'}); % ,{'QV0','Qt','Qx'}
   
   % cone_ndims is the number of space dimensions that the cone actually
   % spans. If the input cone rays define a full-dimenisional cone,
   % cone_ndims is equal to space_ndims. cone_ndims is a natural number of
   % the interval [0, space_ndims].
   cone_ndims=size(conv_cone_k,1);

   if cone_ndims>1 % If cone dimension is at leat 2, we can calculate facets
      % Find all the facets which contain the origin (0). These facets will
      % become the final cone facets
      conv_cone_k_orig_ini=conv_cone_k(:,any(conv_cone_k==0,1));
      n_conv_cone_k_orig=size(conv_cone_k_orig_ini,2);
      if n_conv_cone_k_orig==0 % No facet including the origin
         error('Degenerated input cone? the facets of the input cone cannot be found')
      end
      conv_cone_k_orig_ini(conv_cone_k_orig_ini==0)=[]; % Discard origin point
      % Recover matrix shape after removing origin and substract 1 to account
      % for origin vertex deletion
      conv_cone_k_orig_ini=reshape(conv_cone_k_orig_ini,cone_ndims-1,n_conv_cone_k_orig);
      
      % Find out which rays have been used as vertices of the cone convex hull
      conv_cone_v_ind=unique(conv_cone_k_orig_ini); % conv_cone_v_ind is also sorted by unique()
      % If the origin vertex has been included in the convex hull it should be
      % the first vertex of conv_cone_v_ind, we will this vertex in the final
      % cone vertex matrix
      conv_cone_v=cone_rays_uni(:,conv_cone_v_ind); % Final cone vertex matrix not including the origin
      
      % Translate the facet vertex indices to account for input rays
      % discarded in the resultant convex hull
      conv_cone_k_orig=conv_cone_k_orig_ini;
      removed_cone_v_ind=setdiff(1:size(cone_rays_uni,2),conv_cone_v_ind); % Find out which rays (vertices) have been removed
      for rem_vert_ind=removed_cone_v_ind
         trans_vert_ind=find(conv_cone_k_orig_ini>rem_vert_ind);
         conv_cone_k_orig(trans_vert_ind)=conv_cone_k_orig(trans_vert_ind)-1;
      end
      % conv_cone_k_orig_trans is the translated version of conv_cone_k_orig.
      % So, conv_cone_k_trans is the corresponding K matrix for conv_cone_v (V)
      
      % Triangulate the hypothetical facet that would convert the cone into a
      % polytope (and which does not contain the origin). The vertices of
      % the ontained simplex facets vertices will become the rays of each
      % simplex cone into which the cone convex hull is divided.
      % First, we choose the vertex which is shared by the higher number of
      % facets as common vertex (to lower the number of generated facets):
      % Mode returns the vertex with higher number of repetitions
      % Use it as common vertex for all the simplices that will be generated
      common_vertex=mode(reshape(conv_cone_k_orig, 1, (cone_ndims-1)*n_conv_cone_k_orig));
      % Find all the facets which do not contain the origin
      non_comon_facets=conv_cone_k_orig(:,all(conv_cone_k_orig~=common_vertex,1));
      % Perform the fan triangulation of the facet (polytope)
      n_simplices=size(non_comon_facets,2); % Number of simplices that will be generated
      conv_cone_k_non_orig=[common_vertex*ones(1,n_simplices);non_comon_facets];
   else % 1-dimensional cone
      conv_cone_v_ind=conv_cone_k(:,any(conv_cone_k~=0,1)); % In 1D the facets are vertices
      conv_cone_v=cone_rays_uni(:,conv_cone_v_ind); % Cone vertex matrix: one point or empty
      conv_cone_k_orig=[];
      conv_cone_k_non_orig=[];
   end
else % All the input rays are 0
   conv_cone_v=zeros(space_ndims,0); % Preserve the dimensionality so that total_int() can calculate the integral
   conv_cone_k_orig=[];
   conv_cone_k_non_orig=[];
   cone_ndims=0;
end

% -----------------------------------------------------------------
% At this point:
%  conv_cone_v is the final cone vertex matrix.
%  conv_cone_k_orig is the final cone facet matrix.
%  conv_cone_k_non_orig is the matrix of simplicial cones which compose
%     the cone convex hull.
%  cone_ndims is the number of dimensions of space spanned by the cone.
% From these arrays we generate the output array of cells (elem)
% The content of this cone data structure will be:
% elem{1}: vertex matrix
% elem{2}: elements of dim 1 (rays encoded as vertices): each element is an index to a column of elem{1}
% elem{3}: elements of dim 2 (edges)
% ...
% elem{n_dims-2}: elements of dim. n_dims-3 (peaks)
% elem{n_dims-1}: elements of dim. n_dims-2 (ridges)
% elem{n_dims+0}: elements of dim. n_dims-1 (facets)
% elem{n_dims+1}: elements of dim. n_dims (full-dimensional simplex cones)

% Include vertex matrix (coordinates) in the polytope data structure
elem{1}=conv_cone_v;

% Include vertex elements (also as indices to a vertex matrix column (for homogeneity whic other-dim. elements))
nvertices=size(conv_cone_v,2);
for nvertex=1:nvertices
   elem{2}(nvertex).vertices=nvertex;
   if cone_ndims==2
      elem{2}(nvertex).facets=nvertex; % the next loop will not be executed in 2D, so fill the field now
   else
      elem{2}(nvertex).facets=[]; % The "facets" field will be filled during the next loop or does not need to be filled (1D)
   end
end

% Include facet elements
nfacets=size(conv_cone_k_orig,2);
if cone_ndims>2 % If polytope dimension is lower than 2, facets have already been included (they are vertices)
   for nfacet=1:nfacets
      elem{cone_ndims}(nfacet).vertices=conv_cone_k_orig(:,nfacet); % Facet vertices are not sorted to preserve orientation information
      elem{cone_ndims}(nfacet).facets=nfacet; % Each facet points to itself
      % The index to facet nfacet is added to the corresponding vertex elements (K(:,nfacet))
      old_elem_facets={elem{2}(conv_cone_k_orig(:,nfacet)).facets}; % Previous-iteration facet indices of the involved vertex element (K(:,nfacet))
      new_elem_facets=cellfun(@(el)([el;nfacet]), old_elem_facets, 'UniformOutput', false); % Facet indices of involved vertex element (adding the current facet)
      [elem{2}(conv_cone_k_orig(:,nfacet)).facets]=new_elem_facets{:}; % Update the vertex-element indices to the nfacet facets
   end
end

% Include the simplex cone elements
% The cone decomposition into simplices is defined by conv_cone_k_non_orig.
if cone_ndims>1 % If polytope dimension is lower than 2, facets have already been included (they are vertices)
   n_K_simplices=0; % Number of simplices included
   for n_simplex=1:size(conv_cone_k_non_orig,2)
      % Check that the generated simplices are full dimensional
      if rank(conv_cone_v(:,conv_cone_k_non_orig(:,n_simplex))) == cone_ndims % The simplex has volume, so include it in elem
         n_K_simplices=n_K_simplices+1;
         elem{cone_ndims+1}(n_K_simplices).vertices=conv_cone_k_non_orig(:,n_simplex); % vertices of the simplex cone
         elem{cone_ndims+1}(n_K_simplices).facets=[]; % The cone is not included in any facet
      end
   end
end

% Include the other cone elements
for dim_elem=(cone_ndims-1):-1:3 % For each element type except cone (dim n), facet (dim n-1) and vertex (dim 0): ridge (dim n-2), peak (dim n-3),... edge (dim 1)
   n_faces=[elem{dim_elem+1}.vertices];
   nfaces=size(n_faces,2);
   % Find the facets which share num_shared_vertices vertices
   num_shared_vertices=(dim_elem-2)+1;
   nelem=0; % Number of elements of dim. dim_elem that are being generated
   for nfacet1=1:nfaces
      for nfacet2=(nfacet1+1):nfaces
         shared_vertices=intersect(n_faces(:,nfacet1),n_faces(:,nfacet2));
         if length(shared_vertices) == num_shared_vertices % found at least one element of dimension dim_elem
            % Since the cone must be full-dimensional the number of shared
            % vertices should not be greater than the number of vertices
            % of the searched element.
            
            % Check that that element has not been previously stored
            already_stored=0;
            for nelem2=1:nelem
               if isequal(sort(elem{dim_elem}(nelem2).vertices),sort(shared_vertices))
                  already_stored=nelem2; % Save where the element has been found
                  break
               end
            end
            
            superelement_indices=[nfacet1;nfacet2]; % Indices of the superelements which contain the current element
            if dim_elem < cone_ndims-1 % Translate indices to superelement into indices to facets
               superelements_cells={elem{dim_elem+1}(superelement_indices).facets};  % Indices of the facets which contain the current element
               % Traspose and combine the cell elements to get a column vector
               superelements_cells_t=cellfun(@transpose,superelements_cells,'UniformOutput',false);
               facet_indices=[superelements_cells_t{:}]';
            else % In case the current element if a ridge, we do not need to translate indices
               facet_indices=superelement_indices;
            end
            
            if already_stored==0
               nelem=nelem+1; % Store the element in a new position
               elem{dim_elem}(nelem).vertices=shared_vertices;
               % Store the facet indices in the new element
               elem{dim_elem}(nelem).facets=unique(facet_indices);
            else % Update "facets" field
               elem{dim_elem}(already_stored).facets=unique([elem{dim_elem}(already_stored).facets ; facet_indices]);
            end
         end
      end
   end
end

