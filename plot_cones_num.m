%PLOT_CONES_NUM plots a specified cone and corresponding hypercube 
%   illustrating the squared distance to the cone through a color code.
%   PLOT_CONES_NUM(C,n_pts_dim) plots the cone specified by C and
%   the averaged squared 2-norm residual (norm(d-C*X)^2) for every positive
%   real value of d in the space [0,1]^dims of real numbers. dims is
%   the number of rows of C. X represents the vector that minimizes
%   NORM(d-C*X) subject to X >= 0, which could be calculated by
%   X = LSQNONNEG(C,d). C must be positive and real.
%   The number of columns of C, num_rays, represents the number of vectors
%   passing through 0 (rays) in a dims-dimensional space.
%   [dims,num_rays]=size(C).
%   n_pts_dim specifies the number of points per dimension of the
%   hypercube that will be evaluated and plotted.
%
%   See also LSQNONNEG.

%   Richard R. Carrillo Sánchez. 
%   $Revision: 1.1 $  $Date: 2016/25/01 20:19:00 $
function plot_cones_num(C,n_pts_dim)
plot_flag=1; % Set to 0,1,2 or 3 to control the plot details
minuseps = -eps(100);
plot_red_factor=0.99;

n_dims=size(C,1);
plot_zero_coords=ones(n_dims,1)*(1-plot_red_factor);
center_coords=ones(n_dims,1)*0.5;
C_nonzeros=C(:,any(C)); % remove zero columns
if ~isempty(C_nonzeros)
   % Obtain a full-dimensional cone similar to C to ease calculations
   C_full=fulldim_cone(C_nonzeros);

   % Set the ray length to be confined to the positive orthant, that is,
   % The higher coordinate of every ray must be 1
   C_full=bsxfun(@rdivide,C_full,max(C_full));
   % Discard repeated rays
   cone_v=unique_tol(C_full')'; % The column order may be altered
   
   conv_poly_k=convhulln0(cone_v); % Calculate the facets of the polytope convex hull
    
   if n_dims>1
      % Find all the polytope facets which contain the origin (0).
      % These facets will become the final cone facets
      conv_cone_k=conv_poly_k(:,any(conv_poly_k==0,1));
      % Discard origin vertex
      conv_cone_k(conv_cone_k==0)=[];
      % Recover matrix shape after removing origin
      conv_cone_k=reshape(conv_cone_k,n_dims-1,numel(conv_cone_k)/(n_dims-1));
      % If n_dims is higher than 3, we cannot plot the cone directly.
      % We could calculate the intersection of its extreme rays with an
      % hyperplane, and in case of more than 4 dimensions, project the
      % hyperplane in a 3-dimensional space
      if n_dims>4
         error('Too many dimensions to be shown') % Meanwhile, error
      end
   else
      conv_cone_k=conv_poly_k';
      % Discard origin vertex
      conv_cone_k(conv_cone_k==0)=[];
   end
    
   % Find the extreme rays of the cone
   conv_cone_v=cone_v(:,unique(conv_cone_k));
else % null cone
   cone_v=zeros(n_dims,0);
   conv_cone_v=cone_v;
   conv_cone_k=[];
end

n_rays=size(conv_cone_v,2); % Uses the extreme rays of C as input rays

switch n_dims
   case 0
      hold on
   case 1
      hold on
      for n_ray=1:n_rays
         plot(0,conv_cone_v(1,n_ray),'b.','LineWidth',3)
      end
   case 2
      hold on
      for n_ray=1:n_rays
         plot([0 conv_cone_v(1,n_ray)],[0 conv_cone_v(2,n_ray)],'b','LineWidth',3)
      end
   case 3
      old_font_size=get(0,'DefaultAxesFontSize');
      set(0,'DefaultAxesFontSize',24)
      hold on
      if plot_flag>2
         for n_ray=1:n_rays
            plot3([plot_zero_coords(1) conv_cone_v(1,n_ray)*plot_red_factor],[plot_zero_coords(2) conv_cone_v(2,n_ray)*plot_red_factor],[plot_zero_coords(3) conv_cone_v(3,n_ray)*plot_red_factor],'b','LineWidth',3)
         end
      end
         
      % Draw the cone facets
      for n_facet=1:size(conv_cone_k,2)
         cone_facet=cone_v(:,conv_cone_k(:,n_facet));
            
         edge_intersecs=[];
         % The intersection of the cube edges with the cone facets
         for cube_nedge=1:n_dims
            % Initial vertex of the cube edge
            cube_x=ones(n_dims,1);
            % Direction vector of the current cube edge
            cube_c=zeros(n_dims,1);
            cube_c(cube_nedge)=1;
               
            c=[cone_facet cube_c];
            % Check c for singulary before solving the system
            if rank(c) == min(size(c)) % c has inverse
               % Solve the equation system
               % mu=linsolve(c,cube_x);
               mu=cube_x' / c';
               % Validate the solution (intersection) by checking that:
               % 1- the equation system had 1 solution: all(isfinite(mu))
               % 2- the intersection is inside the current cone facet: all(mu(1:cone_ndim) > minuseps)
               % 3- the intersection is inside the hypercube: all(mu((cone_ndim+1):end) > minuseps) && all(mu((cone_ndim+1):end) < 1+(-minuseps))
               if all(isfinite(mu)) && all(mu > minuseps) && all(mu(3) < 1-minuseps)
                  % Save the intersection point
                  edge_intersecs=[edge_intersecs cube_x-mu(3)*cube_c];
               end
            end
         end
         % Remove duplicate intersection points
         edge_intersecs=unique(edge_intersecs','rows')';
         % Order patch vertices so that a solid polygon is drawn
         if size(edge_intersecs,2)>1 % If more than 1 intersection point
            % If second intersection point is closer to the first vertex than first intersection point
            if norm(edge_intersecs(:,1)-cone_facet(:,1)) > norm(edge_intersecs(:,2)-cone_facet(:,1))
               % Swap points
               edge_intersecs=[edge_intersecs(:,2:end) edge_intersecs(:,1)];
            end
         end
         patch_coords=[cone_facet(:,1) edge_intersecs cone_facet(:,2) plot_zero_coords];
         % Push polygon vertices towards cube center to avoid
         % overlapping with cube faces
         patch_coords=push_vertex(patch_coords,center_coords,plot_red_factor);
         % Draw cone facet
         patch(patch_coords(1,:),patch_coords(2,:),patch_coords(3,:),'red')
      end
      xlabel('x'); ylabel('y'); zlabel('z')
      view(3)
      view(120,35)
      hsv_colormap=hsv(64);
      colormap(hsv_colormap(1:56,:))
      title('error: norm(LSQNN residual)^2')
      colorbar
      daspect([1 1 1])
      box on

      % Allocate space to accelerate execution
      facet_x_coords=zeros(4,n_pts_dim^2);
      facet_y_coords=zeros(4,n_pts_dim^2);
      facet_z_coords(:,:,1)=zeros(4,n_pts_dim^2);
      facet_z_coords(:,:,2)=ones(4,n_pts_dim^2);
      % Calculate coordinates of cube facet squares
      for facet_x=1:n_pts_dim
         for facet_y=1:n_pts_dim
            sq_pos=(facet_x-1)*n_pts_dim + facet_y;
            facet_x_coords(1:4,sq_pos)=[facet_x-1; facet_x; facet_x; facet_x-1]/n_pts_dim;
            facet_y_coords(1:4,sq_pos)=[facet_y-1; facet_y-1; facet_y; facet_y]/n_pts_dim;
         end
      end
      % Calculate the error (color) of the squares
      patch_coords={0,0,0}; % Allocate space for the first patch function arguments (patch coordinates)
      patch_side_coords={facet_x_coords,facet_y_coords}; % Coordinates of the hypercube face
      patch_residual=zeros(1,n_pts_dim^2); % Allocate space for the facet residuals
      if plot_flag > 1
         first_side=1; % Plot all cube facets
      else
         first_side=2; % Only plot 3 cube faces
         % Substitute eliminated faces by plain color ones
         cube_face_xy_coords=[-(1-plot_red_factor) 1*plot_red_factor 1*plot_red_factor -(1-plot_red_factor);-(1-plot_red_factor) -(1-plot_red_factor) 1*plot_red_factor 1*plot_red_factor];
         cube_face_z_coords=ones(1,4)*-(1-plot_red_factor);
         face_color=[0.7 0.7 0.7];
         patch(cube_face_xy_coords(1,:),cube_face_xy_coords(2,:),cube_face_z_coords,face_color);
         patch(cube_face_xy_coords(1,:),cube_face_z_coords,cube_face_xy_coords(2,:),face_color);
         patch(cube_face_z_coords,cube_face_xy_coords(1,:),cube_face_xy_coords(2,:),face_color);
      end
      for ndim=1:3
         dim_list=1:3; % List of hypercube dimensions
         side_dims=setdiff(dim_list,ndim); % Calculate the dimensions (arguments) used to plot the face
         for nside=first_side:2 % 2 faces per dimension
            [patch_coords{side_dims}]=patch_side_coords{:}; % Set those arguments
            patch_coords{ndim}=facet_z_coords(:,:,nside);
            for n_sq=1:n_pts_dim^2 % Evaluate the error of each square
               pt_coords=[mean(patch_coords{1}(:,n_sq)); mean(patch_coords{2}(:,n_sq)); mean(patch_coords{3}(:,n_sq))];
               [~,residual]=lsqnonneg(C,pt_coords);
               patch_residual(n_sq)=residual;
            end
            % delete the squares whose residual is 0
            res_zero_ind=patch_residual<-minuseps;
            patch_residual(res_zero_ind)=[];
            patch_coords{1}(:,res_zero_ind)=[];
            patch_coords{2}(:,res_zero_ind)=[];
            patch_coords{3}(:,res_zero_ind)=[];
            patch(patch_coords{:},patch_residual);
         end
      end
      axis([0 1 0 1 0 1])
      set(0,'DefaultAxesFontSize',old_font_size)
   otherwise
      hold on
      c_plane=1; % hyperplane for intersection with noral vector (1,1,1...) and c=1
      t_plane=c_plane./sum(cone_v,1); % t=c/(xl+yl+zl)
      U_in_plane=cone_v.*(t_plane'*ones(1,n_dims))'; % interseccion=(xl,yl,zl)*t

      v_plane=dir_vector(n_dims); % director vectors used
      v_plane=v_plane(1:end-1,:)'; % vectores directores del plano=[u ; v]
      U_in_plane_zero=U_in_plane-c_plane/n_dims; % Origin of the plane coordinate system=(1,1,1)*c_plane/d_dims
      % 2D n coordinate in the plane=proyection of the 3D vector in the plane n director vector
      % =dot(3D coordinates of the vector,n director vector of the plane)
      % =3D_coor_column'*director_vector
      coord_in_plane=(U_in_plane_zero'*v_plane)';

      scatter3(coord_in_plane(1,:),coord_in_plane(2,:),coord_in_plane(3,:));
      trisurf(conv_cone_k',coord_in_plane(1,:),coord_in_plane(2,:),coord_in_plane(3,:));
      view(3)
      grid on
      daspect([1 1 1])
      title('Intersection of the cone with a 3D space')
end
return


% PUSH_VERTEX moves vertices in v_coords towards point c_coords a factor
% 1-mov_factor
function v_coords2=push_vertex(v_coords,c_coords,mov_factor)
v_coords2=bsxfun(@plus,v_coords*mov_factor,c_coords*(1-mov_factor));


% DIR_VECTOR is an alternative to the Gram–Schmidt process to obtain a list
% of orthogonal vectors
% V=dir_vector(n_dims) returns a matrix of orthogonal vector to the
% vector [1,1,...,1] of dimsension n_dims
function V=dir_vector(n_dims)
V=zeros(n_dims);
nrow=1;
for p_num_com=0:fix(log2(n_dims)-1)
    num_com=2^p_num_com;
    com=[ones(1,num_com) -ones(1,num_com)];
    num_reps=fix(n_dims/length(com));
    num_reps_rem=rem(n_dims,length(com));
    for rep=1:num_reps
        pos=1+(rep-1)*length(com);
        V(nrow,pos:pos+length(com)-1)=com;
        nrow=nrow+1;
    end
    if num_reps_rem >= num_com
        com=[ones(1,num_com*2*num_reps) -2*num_reps*ones(1,num_com)];
        V(nrow,1:1+length(com)-1)=com;
        nrow=nrow+1;
    end
end
V(nrow,:)=ones(1,n_dims);
V_cell=num2cell(V,2);
V_norm=cellfun(@norm,V_cell,'UniformOutput',true)*ones(1,n_dims);
V=V./V_norm; % Vectors in V are normalized
