% elem=devel_cube(n) generates all the elements (hypercube, facets,
% ridges, peaks,..., faces, edges and vertices) of a hypercube of dimension n.
% The hypercube will contain neither the 0 vertex nor any element which
% contain this vertex.
% devel_cube returns elem, which is an array of cells with the format:
% In elem{1} contains V, vertex coordinate matrix.
% elem{ndim+2} which contains the vertex index matrix of the
% elements of dimension ndim. The size of this matrix is ndim+1 x nelems_n
function elem=devel_cube(n_dims)
% Array of cells
% The content of the polytope data structure will be:
% elem{1}: vertex matrix
% elem{2}: elements of dim 1 (vertices): each element is an index to a vertex
% ...
% elem{dim. of the polytpe + 1}: elements of dim. of poly_dims-1 (facets)
% elem{dim. of the polytpe + 2}: element of dim. of poly_dims (cone)
elem={};

% Include vertex matrix (coordinates) in the polytope data structure
% Container of all the vertices of the n_dims-dimensional hypercube
elem{1}=zeros(n_dims,2^n_dims-1); % Allocate space. Each column is a hypercube vertex
for n_vert=uint64(1:(2^n_dims-1)) % All possible vertices except origin
   elem{1}(1:n_dims,n_vert)=double(bitget(n_vert,(1:n_dims)'));
end

% Include vertex elements (also as indices to a vertex matrix column (for homogeneity whic other-dimensionality elements))
elem{2}=1:size(elem{1},2);

% Include higher-dimensionality elements
for dim_elem=2:n_dims % For each element type except vertices: edges, faces,....,ridges and facets
   n_subelems=size(elem{dim_elem},2); % Number of elements of dim. dim_elem-1
   elem{dim_elem+1}=[]; % Create the new cell so that new elements can be appended at the end
   for nelem=1:n_subelems % For each element of dim. dim_elem-1
      subelem_v_ind=elem{dim_elem}(:,nelem);
      subelem_v=elem{1}(:,subelem_v_ind);
      subelem_origin=subelem_v(:,1); % Vertex that will be considered origin of the being-created element
      % Calculates the bitwise OR of all the subelement vertex coordinates
      % subelem_v_or will encode the vertex directions which have been already followed
      subelem_v_or=subelem_origin;
      for nvertex=2:(dim_elem-1)
         subelem_v_or=bitor(subelem_v_or,subelem_v(:,nvertex));
      end
      % Calculate coordinates which are still zero for all the subelement vertices
      % (All possible directions for this elemen: Potential coordinates to increment)
      zeros_coords=find(subelem_v_or==0);
      if ~isempty(zeros_coords)
         % All the potential new elements
         % pot_new_vertex_index=double(bitset(uint64(subelem_v_ind(1)),zeros_coords)); % This is valid only in Matlab
         pot_new_vertex_index=double(arrayfun(@(nbit)(bitset(uint64(subelem_v_ind(1)),nbit)),zeros_coords));
         % Create new elmements with increasing vertices to avoid duplicated elements
         new_vertex_index=pot_new_vertex_index(pot_new_vertex_index > subelem_v_ind(end));
         n_new_elems=length(new_vertex_index); % Number of new elements
         if n_new_elems > 0 % If new elements have been found from the current subelement:
            % Create a matrix containing all the new elments
            % (all the subelement vertices are also included in the new elements)
            new_elems=[repmat(subelem_v_ind,1,n_new_elems);new_vertex_index'];
            % Insert the new elements in elem
            elem{dim_elem+1}=[elem{dim_elem+1} new_elems];
         end
      end
   end
end
