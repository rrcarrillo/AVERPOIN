% int_val=cone_int(cone_elem,sub_v) calculates the integral of the squared
% distance function to a subspace over the polytope defined by the
% intersection points of cone_elem.
% The subspace is defined by the rays (vertices) encoded by the columns
% of sub_v.
% The cone polytope is encoded by cell arrays with the format:
% elem{ndim+1}(nface).vertices which contains the vertices indices of the
% element nface of dimension ndim, except elem{1}(nvertex) which contains
% the vertex coordinate matrix.
% The original cone is encoded by the polytope resulting from its projection on the
% hyperplane with director vector (1,1,...1) (See project_cone()).
function cone_int_val=cone_int(cone_elem, sub_v)
% Simplices whose parallelepiped has less volume that this will be discarted
%vol_eps=eps(2e5); % Ignore volumes lower than those created by fulldim_cone() perturbations
vol_eps=eps(2); % Ignore volumes lower than those created by numerical errors
space_ndims=size(cone_elem{1},1); % number of coordinates of the cone vertices (number of space dimensions in which the cone is)

[cone_k,cone_v]=trian_region(cone_elem); % Decompose into simplices

[n_vertices_facet,n_simplices]=size(cone_k); % number of simplicial cones in which the cone is divided
if n_vertices_facet==space_ndims % full-dimensional region polytope
   int_val=0; % Accumulated sum
   for n_simplex=1:n_simplices % for each simplicial cone
      simp_cone_v=cone_v(:,cone_k(:,n_simplex));
      if rank(simp_cone_v)<space_ndims
         simplex_int=0; % non-full-dimensional simplex
      else
         parallelep_volume=abs(det(simp_cone_v)); % Volume of the corresponding parallelepiped without sign 
         simplex_int=0;
         if parallelep_volume > vol_eps % Otherwise integral simple value will be considered 0 (practically always true)
            simplex_v=[zeros(space_ndims,1) simp_cone_v]; % Add the origin to the simplex vertices
            sub_ortho=gschmidt(sub_v); % Calculate the orthonormal basis of the subspace
            x=sub_ortho*sub_ortho'; % square matrix of size space_ndims x space_ndims
            
            for i1=0:space_ndims
               for i2=i1:space_ndims
                  % The symmetric multilinar form corresponding to the distance
                  % funcion will be evaluated for vertices (rays) v1 and v2
                  v1=simplex_v(:,i1+1);
                  v2=simplex_v(:,i2+1);
                  dist_sym_multi=v1'*v2 - sum(sum((v1*v2' + v2*v1').*x))/2;
                  simplex_int=simplex_int + dist_sym_multi;
               end
            end
            simplex_int=simplex_int*parallelep_volume;
         end
      end
      int_val=int_val+simplex_int;
   end
   % we integrate a 2-homogenous polynomial (q=2)
   % The volume of the simplex is the volume of the parallelep times 1/factorial(n)
   % The formula of the combinations without repetition is
   % factorial(n+q)/(factorial(q)*factorial((n+q)-q), thus we multiply by q!/(q+n)! 
   cone_int_val=int_val*factorial(2)/factorial(2+space_ndims);
else % non full-dimensional region polytope: integral=0
   cone_int_val=0;
end