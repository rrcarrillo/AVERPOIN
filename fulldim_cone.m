% full_rays=fulldim_cone(all_input_rays)
% Calculates a full-dimensional cone from the input cone. The obtained cone
% has similar volume and dimensions to the original cone.
% all_input_rays is a n x r matrix which cotains the coordinates of r rays
% in an n-dimensional space.
% full_rays is a n x rf matrix which cotains the coordinates of rf rays. rf
% if at least n.
function full_rays=fulldim_cone(all_input_rays)
% Only consider non-zero-length rays from input
input_rays=all_input_rays(:,any(abs(all_input_rays) > eps));
% Check than all the rays are in the positive orthant
if any(any(input_rays<0))
   error('All input rays must have positive or zero coordinates')
end

n_dims=size(input_rays,1); % Number of space dimensions

% Normalize ray lengths so that parallel rays can be detected
rays_in_sphere=normalize_vecs(input_rays);

% Remove duplicated rays
unique_rays=unique_tol(rays_in_sphere);

if size(unique_rays,2) > 0 % We deal with the case in wich there are no input rays independently

   % In dimension n a full-dimensional convex polytope has n linearly-independent
   % points plus another point which is different from the others (reference).
   % In this case (a cone) the origin vertex (0) is considered the reference,
   % so we do not need to substract any ray (column) before calculating the
   % matrix rank. In this way, rank returns how many points (rays) are
   % linearly independent.
   % The internal representation used by this software (v-representation)
   % does not deal correctly with the region cones that are generated when
   % the input cone is not full-dimensional, so we make if full-dimensional
   indep_rays=rank(unique_rays); % number of linearly independent rays
   if indep_rays < n_dims % If there are not enough inpependent rays, we add extra rays to make it full dimensional
      % To calculate new rays for the cone we use the base vector of the
      % null space of the cone ray matrix. These vectors are added to the
      % original-cone center ray multiplied by a small factor.
      % To ensure that the resultant new rays are inside the positive
      % orthant we choose the null-space base vector or their opposite,
      % the one which results in a ray nearer to the positive orthant
      % center.

      % Calculate an approximate normal cone center ray
      cone_pcenter=normalize_vecs(mean(unique_rays,2));
      % Calculate the normal center ray of the positive orthant
      orthant_center=normalize_vecs(ones(n_dims,1));
      % Vector used to choose between a null-space base vector or its opposite
      preferred_new_rays_dir=orthant_center-cone_pcenter;

      % Calculate an orthonormal basis for the null space of the ray matrix
      null_space_base=null(unique_rays'); % SVD performed for rank() could have been used for null() for effciency purposes
      % null_space_base contains n_dims-indep_rays vectors which are not
      % collinaear (in fact are normal) to unique_rays
      
      % Calculate the cosine of the angle (scalar product) between
      % preferred_new_rays_dir vectors and the null-space-base vectors.
      % If it is negative, the opposite to the null-space-base vector is
      % more similar to preferred_new_rays_dir
      opposite_null_space_base=null_space_base'*preferred_new_rays_dir < 0;
      
      % Choose between null-space-base vectors and their opposite
      % counterparts. For this, multiply these vector by 1 or -1, depending
      % if opposite_null_space_base values are true (1) or false (0)
      null_space_base_pref=bsxfun(@times,null_space_base,1-2*opposite_null_space_base');
            
      % We do not want the new rays to alter the resultant cone volume,
      % so place them very near (epsilon distance) the original cone center
      % In some cases this value could be increased to improve the volume
      % precision calculation.
      % QHull uses an initial joggle of 30000 times the maximum roundoff
      % error for a distance computation (when option {'QJ'} is specified).
      epsilon=eps(2e5)^(1/(n_dims-indep_rays));
      % new_ray=cone_center+null_vector*epsilon
      epsilon_rays=bsxfun(@plus,null_space_base_pref.*epsilon, cone_pcenter);
      % Append the new rays to the cone ray matrix
      full_rays=[unique_rays normalize_vecs(epsilon_rays)];
      fprintf(1,'New rays added to the initially non-full-dim. cone (extra precision error)\n')
      
      % convhulln works if the number of linearly independent points is higher
      % than the space dimension, so we check it
      if rank(full_rays) < n_dims % At this stage full_rays should encode a full-dimensional cone
         warning('Internal error: a full-dimensional cone could not be obtained from the input rays')
      end
   else
      full_rays=unique_rays; % The input rays already encode a full-dimensional cone
   end
   % At this point the input cone has n_dims linearly independent rays
   % (possibly plus other linearly-dependent rays).
else
   full_rays=zeros(n_dims,0); % Preserve the dimensionality so that total_int() can calculate the integral
end