function plot_example()
C=[2 3 1;3 1 1;0 0 1]';

n_dims=size(C,1); % Num. of ray coordinates

cube_elem=devel_cube(n_dims); % Calculates the elements (vertices, faces,...) of an n_dims-dimensional hypercube
full_cone_c=fulldim_cone(C); % Obtain a full-dimensional cone similar to cone_rays (this is only needed for non-full-dimensional input cones)
cone_elem=devel_cone(full_cone_c); % Calculates the elements of full_cone_c
cone_elem=add_regions(cone_elem); % adds a region cone for each element of cone_elem
cone_elem=cone_intersec(cone_elem,cube_elem); % calculates the intersection of cone_elem cone with the hypercube
cone_elem=region_intersec(cone_elem,cube_elem); % calculates the intersection of each region cone with the hypercube

adj_cone_elem=cone_elem{2}(1);
adj_facets=[cone_elem{3}(adj_cone_elem.facets).normal];
adj_cone = [adj_facets cone_elem{1}(:,adj_cone_elem.vertices)];

C=normc(C);
adj_cone=normc(adj_cone);

plot_cube(false)
plot_rays(C, 'black')
pause

plot_cube(false)
plot_rays(C, 'black')
plot_edges(C, 'red')
plot_hull(C, [0 0.7 0.7], false)
pause 

plot_cube(false)
plot_rays(C, 'black')
plot_hull(C, [0 0.7 0.7], false)
plot_rays(adj_cone, 'black')
pause

plot_cube(false)
plot_rays(C, 'black')
plot_hull(C, [0 0.7 0.7], false)
plot_rays(adj_cone, 'black')
plot_edges(adj_cone, 'red')
plot_hull(adj_cone, [1 0.5 0], false)
pause

plot_cube(true)
plot_rays(C, 'black')
plot_hull(C, [0 0.7 0.7], false)
plot_rays(adj_cone, 'black')
plot_hull(adj_cone, [1 0.5 0], false)
[cone_k,cone_v]=trian_region(adj_cone_elem.region);
plot_points(cone_v, 'red')
pause

plot_cube(true)
plot_rays(C, 'black')
plot_hull(C, [0 0.7 0.7], false)
%plot_rays(adj_cone, 'black')
plot_hull(cone_v, [1 0.5 0], true)
%pause

function plot_points(V, col)
plot3(V(1,:),V(2,:),V(3,:),[col 'o'],'LineWidth',2)

function plot_hull(V, col, plot_base)
[space_ndims,~]=size(V);

V_inter=zeros(space_ndims,1); % Include origin point
V_inter=[V_inter V];

K_inter=convhulln0(V_inter(:,2:end))'+1;

if ~plot_base
   K_inter=K_inter(any(K_inter == 1,2),:); % discard non-conical faces
end

trisurf(K_inter,V_inter(1,:)*1,V_inter(2,:)*1,V_inter(3,:)*1, 'Facecolor', col, 'LineWidth',1);

function plot_rays(V, col)
[~,cone_nrays]=size(V);

for ray_ind=1:cone_nrays
   quiver3(0,0,0,V(1,ray_ind)*1.3,V(2,ray_ind)*1.3,V(3,ray_ind)*1.3,col,'LineWidth',2)
end

function plot_edges(V, col)
[~,cone_nrays]=size(V);

for ray_ind=1:cone_nrays
   plot3([0 V(1,ray_ind)],[0 V(2,ray_ind)],[0 V(3,ray_ind)],col,'LineWidth',6)
end

function plot_cube(full_cube)
clf
hold on
view(3)
axis vis3d % Disable stretch-to-fill
%grid on
axis off
daspect([1 1 1])
%axis([-1 1 -1 1 -1 1])
view(160, 32)

V_2face=de2bi([0 1 3 2],2)';

V_3face0=[V_2face; zeros(1,4); V_2face];
V_3face1=[V_2face; ones(1,4); V_2face];

for n_face=1:3
   facet_c = V_3face0(n_face:(n_face+2),:);
   patch(facet_c(1,:),facet_c(2,:),facet_c(3,:),'red','FaceAlpha',0,'EdgeColor',[0.7 0.7 0],'LineWidth',2)
   if full_cube
      facet_c = V_3face1(n_face:(n_face+2),:);
      patch(facet_c(1,:),facet_c(2,:),facet_c(3,:),'blue','FaceAlpha',0,'EdgeColor',[0.7 0.7 0],'LineWidth',2)
   end
end


