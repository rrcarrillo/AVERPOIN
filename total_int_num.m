%TOTAL_INT_NUM Average value of the squared 2-norm of the minimization
%   residual using numerical integration through the rectangle method.
%   res=TOTAL_INT_NUM(C,numerical_resolution_in_points) approximates the
%   averaged squared 2-norm residual (norm(d-C*X)^2) for a set of
%   positive real values of d in the space [0,1]^dims of real numbers.
%   dims is the number of rows of C. X represents the vector that minimizes
%   NORM(d-C*X) subject to X >= 0, which could be caluclated by
%   X = LSQNONNEG(C,d). C must be positive and real.
%
%   The number of columns of C, num_rays, represents the number of vectors
%   passing through 0 (rays) in a dims-dimensional space.
%   [dims,num_rays]=size(C).
%   the averaged residual is calculated by through the squared 2-norm
%   distance from a point v in the dims-dimensional space to the convex
%   hull of the cone defined by the column vectors of C over the [0,1]^dims
%   hypercube.
%   numerical_resolution_in_points specifies the number of points that
%   will be avaluated per dimension.
%   res fields will contain the following approximate values:
%    res.cone_volume: Very approximate volume of the cone intersection with
%                      the hypercube
%    res.integ: Average value of the squared 2-norm of the minimization
%                      residualthe (integration of this fn over the hypercube)
%    res.max_integ: Maximum value that the integration can reach in n dimensions
%    res.norm_integ: Normalized the integration value
%    res.fitness: 1-res.norm_integ (an opposite measurement)
%
%   See also RANK, LSQNONNEG.

%   Richard R. Carrillo Sánchez.
%   $Revision: 1.1 $  $Date: 2017/04/15 16:43:00 $

% References:
% - Lawson and Hanson, "Solving Least Squares Problems", Prentice-Hall, 1974.
% - J. A. De Loera, B. Dutra, M. KöPpe, S. Moreinis, G. Pinto, and J. Wu.
%     2013. Software for exact integration of polynomials over polyhedra.
%     Comput. Geom. Theory Appl. 46, 3 (April 2013), 232-252.
function res=total_int_num(C,num_point_per_dim)
large_eps=sqrt(eps(1));
C_nonzeros=C(:,any(C)); % remove zero columns
% Prevent C from being empty
if isempty(C_nonzeros)
    C_nonzeros=zeros(size(C,1),0);
else
    % Set the ray length to be confined to the positive orthant, that is,
    % the higher coordinate of every ray must be 1
    C_nonzeros=bsxfun(@rdivide,C_nonzeros,max(C_nonzeros));
    % Discard repeated rays
    C_nonzeros=unique(C_nonzeros','rows')'; % The column order may be altered
end

n_dims=size(C_nonzeros,1); % Uses the shrinked C as algorithm input

max_coord_value=1; % maximum value of a coordinate
coord_increment=max_coord_value/num_point_per_dim; % the coordinate is incremented this amount (h) each iteration
residual_num_integ=0; % Numerical integration (sum) of the squared residuals
cone_vol_num_integ=0; % Numerical integration (sum) of the cone volume (inside the hypercube)
int_cell_volume=coord_increment^n_dims; % Volume of the integration cell (rectangle)
%total_cube_vol=0; % Total hypercube volume (for calculation checking)

elem_coords=ones(n_dims,1)*coord_increment/2; % current integration-point coordinate column vector: the center of the rectangle top is set by the the int. fn.
n_cur_point=1; % number of the current point
increment_end=false;
while increment_end==false
    [~,residual]=lsqnonneg(C,elem_coords);
    residual_num_integ=residual_num_integ+residual*int_cell_volume;
    cone_vol_num_integ=cone_vol_num_integ+(residual<large_eps)*int_cell_volume;
    %total_cube_vol=total_cube_vol+int_cell_volume;
    
    [elem_coords,increment_end]=inc_coords(elem_coords,max_coord_value,coord_increment);
    n_cur_point=n_cur_point+1;
end

res.cone_volume=cone_vol_num_integ; % Very approximate
%res.volume_error=1-total_cube_vol;
res.max_integ=n_dims/3;
res.integ=residual_num_integ;
res.norm_integ=res.integ/res.max_integ;
res.fitness=1-res.norm_integ;

% INC_COORDS Increments the Cartesian coordinates of a point over
% all the space
function [coords,increment_end]=inc_coords(coords,max_coord_value,coord_increment)
ndims=length(coords);
increment_end=true;
for ncompon=1:ndims
    if coords(ncompon)+coord_increment < max_coord_value % the max coord value may not actualy be evaluated
        coords(ncompon)=coords(ncompon)+coord_increment;
        increment_end=false;
        break;
    else
        coords(ncompon)=coord_increment/2; % reset coordinate to h/2 (rectangle height is calculated at the center of the top face: Midpoint approximation)
    end
end
