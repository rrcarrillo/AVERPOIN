%INT_RES Average value of the squared 2-norm of the minimization residual.
%   res=int_res(C,plot_flag) calculates the averaged squared 2-norm
%   residual NORM(d-C*X)^2 for every positive real value of d in the
%   space [0,1]^n of real numbers. n is the number of rows of C.
%   X represents the vector that minimizes NORM(d-C*X) subject to X >= 0,
%   which could be caluclated by X = LSQNONNEG(C,d). C must be positive
%   and real.
%   plot_flag is an optional argument which indicates if the n-dimenional
%   hypercube and cone (or its projection if n>4) must be plotted (true) to
%   illustrate the aspect of the convex hull.
%   res is the output structre which will contain the values calculated by
%   INT_RES in the following fields:
%    res.volume_error: Numeric error caused by the division of the
%                      hypercube volume into integration regions.
%    res.cone_volume: Volume of the cone intersection with the hypercube
%    res.integ: Average value of the squared 2-norm of the minimization
%                      residualthe (integration of this fn over the
%                      hypercube). Ir().
%    res.max_integ: Maximum value that the integration can reach in n
%                      dimensions.
%    res.norm_integ: Normalized the integration value. IrN().
%    res.fitness: 1-res.norm_integ; % an opposite measurement. 1-IrN().
%
%   The number of columns of C, r, represents the number of vectors
%   from 0 (rays) in a n-dimensional space. [n,r]=size(C).
%   This averaged residual is calculated by integrating the squared 2-norm
%   distance from a point v in the dims-dimensional space to the convex
%   hull of the cone defined by the column vectors of C over the [0,1]^n
%   hypercube.
%
%   See also RANK, LSQNONNEG.

%   This file is part of AVERPOIN
%   (AVErage of squared Residuals for Positive-Input Networks)
%   $Revision: 1.5 $  $Date: 2015/07/12 20:19:00 $
%     Copyright (C) 2015  Richard R. Carrillo (University of Granada)
% 
%     AVERPOIN is free software: you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as
%     published by the Free Software Foundation, either version 3 of the
%     License, or (at your option) any later version.
% 
%     AVERPOIN is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%     GNU Lesser General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program. If not, see <http://www.gnu.org/licenses/>.

% References:
% - Lawson and Hanson, "Solving Least Squares Problems", Prentice-Hall, 1974.
% - J.B. Lasserre, K.E. Avrachenkov, The multi-dimensional version of
%   int(a,b,x^p)dx, Amer. Math. Mon. 108 (2001) 151–154.
function res=int_res(cone_rays,plot_flag)
n_dims=size(cone_rays,1); % Num. of ray coordinates

if n_dims >= 8
   warning('This program has been developed mainly for demonstration purposes and it is not optimized. Use a cone_rays input matrix of less than 8 rows to obtain the result in a reasonable time.')
end

cube_elem=devel_cube(n_dims); % Calculates the elements (vertices, faces,...) of an n_dims-dimensional hypercube

full_cone_c=fulldim_cone(cone_rays); % Obtain a full-dimensional cone similar to cone_rays (this is only needed for non-full-dimensional input cones)
cone_elem=devel_cone(full_cone_c); % Calculates the elements of full_cone_c
cone_elem=add_regions(cone_elem); % adds a region cone for each element of cone_elem

cone_elem=cone_intersec(cone_elem,cube_elem); % calculates the intersection of cone_elem cone with the hypercube
cone_elem=region_intersec(cone_elem,cube_elem); % calculates the intersection of each region cone with the hypercube

vol_list=total_vol(cone_elem); % Volume of region simplicies: the first element is the volume of the cone simplices
res.cone_volume=vol_list(1); % volume of the cone simplices
res.volume_error=1-sum(vol_list); % sum the volume of each simplex into which the hypercube has been divided
res.max_integ=n_dims/3; % Maximum value that the integral can reach in n_dims dimensions
res.integ=total_int(cone_elem); % sum the integral value of the squared distance over each region simplex
res.norm_integ=res.integ/res.max_integ; % normalize the integral value
res.fitness=1-res.norm_integ; % measurement opposite to the integral value (fitness)

if abs(res.volume_error) > 1e-5
   warning(['The calculated integral value could be imprecise. Vol. error: ' num2str(res.volume_error)])
end

if nargin>1 && plot_flag
   plot_cones(cone_elem) % Plot the regions in the hypercube
end