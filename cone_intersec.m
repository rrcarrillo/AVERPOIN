% cone_elem=cone_intersec(cone_elem,cube_elem) calculates the intersections
% of all the cone elements (cone, facets, ridges, peaks,..., faces, edges and
% vertices) with their correponding hypercube element (vertices, edges,..., facets).
% The two polytopes are encoded by cell arrays with the format:
% elem{ndim+1}(nface).vertices which contains the vertices indices of the
% element nface of dimension ndim, except elem{1} which contains the
% vertex coordinate matrix.
% The cone is encoded by the polytope resulting from its projection on the
% hyperplane with director vector (1,1,...1) (See project_cone()).
function cone_elem=cone_intersec(cone_elem,cube_elem)
minuseps = -eps(100);
cone_ndims=length(cone_elem)-1; % Dimensions of the space that the cone spans
cube_ndims=length(cube_elem)-1; % Dimensions of the space that the hypercube spans
if cone_ndims > cube_ndims
   warning('The number of dimensions of the space spanned by the cone should be the same or lower lower than the hypercube ones')
end
cone_vertices=cone_elem{1}; % Vertices (rays) of the cone (which in the cone are considered dim. 1)
cube_vertices=cube_elem{1}; % Vertices of the hypercube (which in the cone are considered dim. 0)
space_ndims=size(cube_vertices,1); % Check the number of coordinates of the first cone vertex to find out the space number of dimensions

for cone_ndim=1:cone_ndims % For each type of cone element: rays, faces,..., cone
   % In the cone, cone_elem{2} are vertices (rays) (dim. 1)
   cube_ndim=space_ndims-cone_ndim; % The elements of the hypercube (cone_elem{cube_ndim+1}) which will be checked for intersection have dim. cube_ndim (cube_elem{cube_ndim+2})
   cone_nelems=length(cone_elem{cone_ndim+1}); % Num. of cone elems. of this dim.
   for cone_nelem=1:cone_nelems % For each element of this type
      if ~isfield(cone_elem{cone_ndim+1}(cone_nelem),'intersections') % Create the filed if it didn't exist
         cone_elem{cone_ndim+1}(cone_nelem).intersections=[];
      end
      % Find out the (vertices) rays of the current cone element
      cone_c=cone_vertices(:,cone_elem{cone_ndim+1}(cone_nelem).vertices)';
      actual_cone_elem_dim=size(cone_c,1);
      if cone_ndim > actual_cone_elem_dim % The current element has lower dimnesion than it should have
         continue % Skip the current element
      end

      cube_nelems=size(cube_elem{cube_ndim+2},2); % Num. of cube elems. for this dim.
      % The intersection of the hypercube vertices with the cone
      % space is considered apart, so, hypercube vertices are
      % not considered alone for intersection in this loop
      for cube_nelem=1:cube_nelems
         % Vertices of the current hypercube element
         cube_facet=cube_vertices(:,cube_elem{cube_ndim+2}(:,cube_nelem));
         % Define the hypercube facet as: cube_x+mu_1*v_2+...+mu_n*v_(n+1).
         % The hypercube n-facet has 2^n vertices, we only need n+1 (cube_ndim)
         % vertices: discard the rest.
         [~,cube_facet_order]=sort(sum(cube_facet~=0,1)); % The vertices with less non-zero coordinates are choosen first
         % Assume that cube_x is the vertex nearest to the hypercube origin
         cube_x=cube_facet(:,cube_facet_order(1))';
         % Substract cube_x (hypercube facet origin) from the other used facet vertices
         cube_c=bsxfun(@minus,cube_facet(:,cube_facet_order(2:(cube_ndim+1)))',cube_x);
         % The choosen vertices in cube_c must be perpendicular in order to
         % simplify the hit test in the hypercube facet
         cube_c2=cube_c*cube_c'; % Check that only diagonal elements are different from 0
         cube_c2(logical(eye(length(cube_c2))))=0;
         if ~all(all(cube_c2==0))
            warning('Internal error: Non-perpendicular vertors: Modify algorithm for selecting hypercube facet vertices')
         end
         c=[cone_c;-cube_c];
         % Check c for singulary before solving the system
         % We must avoid using det() for this when using floating point computations
         if rank(c) == min(size(c)) % c has inverse
            % Solve the equation system
            mu=cube_x / c;
            % Validate the intersection by checking that:
            % 1- the equation system had 1 solution: all(isfinite(mu))
            % 2- the intersection is inside the current cone facet: all(mu(1:cone_ndim) > minuseps)
            % 3- the intersection is inside the hypercube: all(mu((cone_ndim+1):end) > minuseps) && all(mu((cone_ndim+1):end) < 1+(-minuseps))
            if all(isfinite(mu)) && all(mu > minuseps) && all(mu((actual_cone_elem_dim+1):end) < 1-minuseps)
               % Save the intersection point
               inter_point=mu(1:actual_cone_elem_dim)*cone_c;
               cone_elem{cone_ndim+1}(cone_nelem).intersections=[cone_elem{cone_ndim+1}(cone_nelem).intersections inter_point'];
            end
         end
      end
      % Remove duplicate intersection points
      cone_elem{cone_ndim+1}(cone_nelem).intersections=unique_tol(cone_elem{cone_ndim+1}(cone_nelem).intersections);
   end   
end
