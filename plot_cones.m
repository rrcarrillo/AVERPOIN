% plot_cones(cone_elem)
% plots the cone defined by cone_elem through cell arrays with the format:
% elem{ndim+1}(nface).vertices which contains the vertices indices of the
% element nface of dimension ndim, except elem{1} which contains the
% vertex coordinate matrix.
% The cone is encoded by the polytope resulting from its projection on the
% hyperplane with director vector (1,1,...1) (See project_cone()).

% Copyright (C) 2018 Richard R. Carrillo (University of Granada)
% This file is part of AVERPOIN. See int_res.m.
function plot_cones(cone_elem)
V=cone_elem{1};
[space_ndims,cone_nrays]=size(V); % Number of space dimensions and number of rays
if space_ndims>=2
   figure
   set(gcf,'Name',mat2str(V,2))
   set(gca, 'FontSize', 14)
end
switch space_ndims
   case 0 % 0D cone
   case 1 % 1D cone
   case 2 % 2D cone
      cone_ndims=length(cone_elem)-2;

      box on
      axis([0 1 0 1])
      xlabel('x')
      ylabel('y')
      daspect([1 1 1])
      grid on
      title('2D cone')
      hold on
      
      % plot_cone(cone_elem, 5); % Plot initial cone
      % plot rays
      for ray_ind=1:cone_nrays
         plot([0 V(1,ray_ind)],[0 V(2,ray_ind)],'k','LineWidth',3)
      end
      drawnow

      % Plot region cones
      ncolor=0;
      for cone_ndim=1:cone_ndims % For each type of cone element: rays, faces,..., except cone
         cone_nelems=length(cone_elem{cone_ndim+1}); % Num. of cone elems. of this dim.
         for cone_nelem=1:cone_nelems % For each element of this type
            plot_cone(cone_elem{cone_ndim+1}(cone_nelem).region, mod(ncolor,7)+1);
            ncolor=ncolor+1;
         end
      end

      hold off
   case 3 % 3D space
      cone_ndims=length(cone_elem)-2;

      box on
      axis([-0.1 1 -0.1 1 -0.1 1])
      xlabel('x')
      ylabel('y')
      zlabel('z')
      daspect([1 1 1])
      view(3)
      grid on
      axis vis3d
      title('3D cone')
      view(150,36)
      hold on
      
      % plot_cone(cone_elem, 5); % Plot initial cone
      
      % Plot initial-cone rays
      for ray_ind=1:cone_nrays
         plot3([0 V(1,ray_ind)],[0 V(2,ray_ind)],[0 V(3,ray_ind)],'c','LineWidth',3)
      end
      drawnow

      % Plot region cones
      ncolor=0;
      for cone_ndim=1:cone_ndims % For each type of cone element: rays, faces,..., except cone
         cone_nelems=length(cone_elem{cone_ndim+1}); % Num. of cone elems. of this dim.
         for cone_nelem=1:cone_nelems % For each element of this type
            plot_cone(cone_elem{cone_ndim+1}(cone_nelem).region, mod(ncolor,7)+1);
            ncolor=ncolor+1;
         end
      end

      hold off
   case 4 % 4D space
      cone_ndims=length(cone_elem)-2;
      
      hold on
      view(3)
      axis vis3d % Disable stretch-to-fill
      grid on
      daspect([1 1 1])
      title('Intersection of the 4D cone with the plane with normal vector (1,1,1)')

      % Plot cone
      coord_in_hplane=project_on_hyperplane(V,3);
      K_in_plane=[cone_elem{end-1}.vertices];
      trisurf(K_in_plane',coord_in_hplane(1,:),coord_in_hplane(2,:),coord_in_hplane(3,:), 'Facecolor','k', 'EdgeColor','w');
      text(coord_in_hplane(1,:),coord_in_hplane(2,:),coord_in_hplane(3,:),num2str([1:cone_nrays]'), 'BackgroundColor','b', 'Color', 'w');
      drawnow

      % Plot region cones
      ncolor=0;
      for cone_ndim=1:cone_ndims % For each type of cone element: rays, faces,..., except cone
         cone_nelems=length(cone_elem{cone_ndim+1}); % Num. of cone elems. of this dim.
         for cone_nelem=1:cone_nelems % For each element of this type
            plot_cone(cone_elem{cone_ndim+1}(cone_nelem).region, mod(ncolor,7)+1);
            ncolor=ncolor+1;
         end
         drawnow
      end
      
      hold off
   otherwise % if cone ndims > 4, project in 3D space
      cone_ndims=length(cone_elem)-2;
      
      hold on
      view(3)
      axis vis3d % Disable stretch-to-fill
      grid on
      daspect([1 1 1])
      title('Intersection of the cone with the (1,1,...1) hyperplane and projection on 3D subspace')

      % Plot cone
      coord_in_hplane=project_on_hyperplane(V,3);
      K_in_plane=[cone_elem{end-1}.vertices];
      trisurf(K_in_plane',coord_in_hplane(1,:),coord_in_hplane(2,:),coord_in_hplane(3,:), 'Facecolor','k', 'EdgeColor','w');
      text(coord_in_hplane(1,:),coord_in_hplane(2,:),coord_in_hplane(3,:),num2str([1:cone_nrays]'), 'BackgroundColor','b', 'Color', 'w');
      drawnow

      % Plot region cones
      ncolor=0;
      for cone_ndim=1:cone_ndims % For each type of cone element: rays, faces,..., except cone
         cone_nelems=length(cone_elem{cone_ndim+1}); % Num. of cone elems. of this dim.
         for cone_nelem=1:cone_nelems % For each element of this type
            plot_cone(cone_elem{cone_ndim+1}(cone_nelem).region, mod(ncolor,7)+1);
            ncolor=ncolor+1;
         end
         drawnow
      end
      
      hold off
end

% plot_cone(cone_elem, face_color)
function plot_cone(cone_elem, face_color)
V=cone_elem{1};
[space_ndims,cone_nrays]=size(V); % Number of space dimensions and number of rays

color_names={'yellow','magenta','cyan','red','green','blue','white','black'};

cone_ndims=length(cone_elem)-2; % Dimensions of the space that the (projected) cone spans
V_inter=zeros(space_ndims,1); % Include origin point
for cone_ndim=1:(cone_ndims+1) % For each type of cone element: rays, faces,..., cone
   V_inter=[V_inter cone_elem{cone_ndim+1}.intersections]; % Collect all the intersection points
end
V_inter=unique_tol(V_inter); % Non-full-dimensional polytopes will not be plotted

hold on

switch space_ndims
   case 2
      K_inter=convhulln0(V_inter(:,2:end))'+1;
      K_inter_lin=reshape(K_inter',1,numel(K_inter));
      % * Algorithm to arrange the points so they are in counterclockwise
      %   which is needed to draw the polygon (MATLAB convhulln does it)
      V_conv=V_inter(:,K_inter_lin);
      % Calculate an approximate polygon center ray
      V_conv_center=mean(V_conv,2);
      % Translate points to the approx. polygon center
      V_conv_trans=bsxfun(@minus,V_conv,V_conv_center);
      % Calculate point angles to the x axis
      V_conv_angle=atan2(V_conv_trans(2,:),V_conv_trans(1,:));
      % Sort angles counterclockwise
      [~,V_sort_ind]=unique(V_conv_angle);
      % Repeat last point
      V_sort=[V_conv(:,V_sort_ind) V_conv(:,V_sort_ind(1))];
      % Draw polygon
      patch(V_sort(1,:),V_sort(2,:),color_names{face_color})
      %for ray_ind=1:cone_nrays
      %   plot([0 V(1,ray_ind)],[0 V(2,ray_ind)],color_names{face_color},'LineWidth',2)
      %end
   case 3
      K_inter=convhulln0(V_inter(:,2:end))'+1;
      %K_inter=convhullan(V_inter(:,1:end)'); % convhullan() can be used to plot non-full-dimensional polytopes
      trisurf(K_inter,V_inter(1,:),V_inter(2,:),V_inter(3,:), 'Facecolor',color_names{face_color});
      %for ray_ind=1:cone_nrays
      %   plot3([0 V(1,ray_ind)],[0 V(2,ray_ind)],[0 V(3,ray_ind)],color_names{face_color},'LineWidth',3)
      %end
   case 4
      V_inter_hplane=project_on_hyperplane(V_inter(:,2:end),3);
      K_inter=convhulln0(V_inter(:,2:end)); % Skip the origin point (it is included by convhulln0())
      
      % Find all the facets which contain the origin (0). These facets will
      % become the finally plotted facets
      K_inter_ini=K_inter(:,any(K_inter==0,1));
      K_inter_size=size(K_inter_ini);
      K_inter_ini(K_inter_ini==0)=[]; % Discard origin point
      K_inter_ini=reshape(K_inter_ini,K_inter_size(1)-1,K_inter_size(2)); % Recover matrix shape

      trisurf(K_inter_ini',V_inter_hplane(1,:),V_inter_hplane(2,:),V_inter_hplane(3,:), 'Facecolor',color_names{face_color}, 'facealpha',0.2);
   otherwise
      V_inter_hplane=project_on_hyperplane(V_inter(:,2:end),3); % intersect and proyect
      V_inter_hplane_trans=bsxfun(@minus,V_inter_hplane(:,2:end),V_inter_hplane(:,1)); % Translate to the origin (convhulln0 only works with polytopes with this vertex)
      K_inter=convhulln0(V_inter_hplane_trans)'+1;
      % Due to the region proyections on a lower dimensional space, plotted
      % polyhedra will probably overlap
      trisurf(K_inter,V_inter_hplane(1,:),V_inter_hplane(2,:),V_inter_hplane(3,:), 'Facecolor',color_names{face_color}, 'facealpha',0.1);
end

% Project the vertices defined by the columns of V on the hyperplane whose
% normal vector is (1,1,....,1) and spans a space of ndims_hyper
function coord_in_plane=project_on_hyperplane(V,ndims_hyper)
space_ndims=size(V,1); % Number of space dimensions
% Intersection of the hyperplane with noral vector (1,1,1...) and c=1 with
% the column vectors of cone_rays (u vectors)
c_plane=1;
t_plane=c_plane./sum(V,1); % t=c/(xl+yl+zl)
U_in_plane=V.*(t_plane'*ones(1,space_ndims))'; % interseccion=(xl,yl,zl)*t
v_plane=dir_vector(space_ndims);
v_plane=v_plane(:,1:end-(space_ndims-ndims_hyper)); % vectores directores del plano=[u ; v]
U_in_plane_zero=U_in_plane-c_plane/space_ndims; % Origin of the plane coordinate system=(1,1,1)*c_plane/d_dims
% 2D n coordinate in the plane=proyection of the 3D vector in the plane n director vector
% =dot(3D coordinates of the vector,n director vector of the plane)
% =3D_coor_column'*director_vector
coord_in_plane=(U_in_plane_zero'*v_plane)';

% A fast alternative to the Gram�Schmidt process to obtain a list of
% orthogonal vectors
% V=dir_vector(n_dims) returns a matrix of column vectors which are
% orthogonal to the vector (1,1,...,1) of dimsension n_dims. This vector is
% also included in the output matrix, so the size of V is n_dims x n_dims.
% The columns of V are normalized so that, V become an orthonormal basis
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
V=(V./V_norm)'; % Vectors in V are normalized
