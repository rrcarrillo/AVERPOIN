function cone_elem=prureg(cone_c)
%cone_c=[0.8 1 0.2; 0.2 1 0.2; 0.5 1 0.8]'; % triangulo 3D en cara
%cone_c=[1 0.5 0.5; 0.5 1 0.5; 0.5 0.5 1]'; % triangulo 3D en esquina

%cone_c=[0.8 1 0.2; 0.2 1 0.2; 0.2 1 0.8; 0.8 1 0.8]'; % cuadrado 3D
%cone_c=[1 0.8 0.3; 0 0.8 0.2; 0.2 1 0.8; 1 0.8 1]'; % cuadrado 3D grande
%cone_c=[0.2 1 0.2 0.2; 0.5 0.9 1 0.2; 0.8 1 0.2 0.2; 0.5 1 0.5 0.8]'; % 4D
%cone_c=[1 1 0]';
%cone_c=[0.2 1 0.2 1; 0.8 1 0.2 1; 0.2 1 0.8 1; 0.8 1 0.8 1]'; % degen
%cone_c=rand(5,6);

%cone_c=[1 0.2 0.2 0.2 0.2; 1 0.8 0.2 0.2 0.2; 1 0.2 0.8 0.2 0.2; 1 0.8 0.8 0.2 0.2; 1 0.2 0.2 0.8 0.2; 1 0.2 0.2 0.2 0.8]'; % ; cuadrado 5D full-dimensional con dos facetas coplanares
%cone_c=[1 0.2 0.2 0.2; 1 0.8 0.2 0.2; 1 0.2 0.8 0.2; 1 0.8 0.8 0.2; 1 0.2 0.2 0.8]'; % cuadrado 4D full-dimensional con dos facetas coplanares

%cone_c=[0 1; 1 0]';
%cone_c=[1 1]';

%cone_c=[0.8 1 0.2; 0.2 1 0.2; 0.5 1 0.3]';

% cone_c=[0.1  0.2  0.3  0.4  0.5;
%         0.9  0.5  0.1  0.8  0.7;
%         0.8  0.3  0.3  0.5  0.3;
%         0.7  0.8  0.5  0.6  0.5;
%         0.2  0.3  0.7  0.2  0.8;
%         0.1  0.2  0.9  0.2  0.4]'; % wrong volume calc.
% cone_c=[0.1294 0.2757 0.3868 0.4359 0.5468;
%         0.9063 0.5085 0.5200 0.8176 0.7948;
%         0.8443 0.3786 0.3116 0.5328 0.3507;
%         0.7390 0.8759 0.5502 0.6225 0.5870;
%         0.2077 0.3012 0.7709 0.2305 0.8443;
%         0.1948 0.2259 0.9707 0.2277 0.4357]'; % Wrong volume calc.

n=size(cone_c,1);

cube_elem=devel_cube(n);

full_cone_c=fulldim_cone(cone_c);
cone_elem=devel_cone(full_cone_c);
cone_elem=add_regions(cone_elem);

cone_elem=cone_intersec(cone_elem,cube_elem);
cone_elem=region_intersec(cone_elem,cube_elem);

total_vol_comp=1-sum(total_vol(cone_elem));
max_int=n/3;
total_int_val=total_int(cone_elem);

cone_whos=whos('cone_elem');
fprintf(1,'Integ.: %g of %g (Vol. error: %g)\n',total_int_val, max_int, total_vol_comp)
fprintf(1,'Normal. integ.: %g (fitness: %g)\n',total_int_val/max_int, 1-total_int_val/max_int)
%fprintf(1,'Cone data struct. size: %ibytes\n',cone_whos.bytes)
plot_cones(cone_elem)
