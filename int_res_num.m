%INT_RES_NUM Numerical version of the INT_RES function for testing and
% illustration. It approximates numerically the average value of the
%   squared 2-norm of the minimization residual.
%   res=int_res_num(C,n_pts_dim,plot_flag) calculates the averaged squared
%   2-norm residual NORM(d-C*X)^2 for n_pts_dim^n real values of d in the
%   space [0,1]^n of real numbers. n is the number of rows of C and
%   n_pts_dim specifies how many values will be avaluated per dimension.
%   X represents the vector that minimizes NORM(d-C*X) subject to X >= 0,
%   which could be caluclated by X = LSQNONNEG(C,d). C must be positive
%   and real.
%   plot_flag is an optional argument which indicates if the n-dimenional
%   hypercube and cone (or its projection if n>4) must be plotted (true) to
%   illustrate the aspect of the convex hull.
%   res fields will contain the following approximate values:
%    res.cone_volume: Volume of the cone intersection with the hypercube
%    res.integ: Average value of the squared 2-norm of the minimization
%                      residualthe (integration of this fn over the hypercube)
%    res.max_integ: Maximum value that the integration can reach in n dimensions
%    res.norm_integ: Normalized the integration value
%    res.fitness: 1-res.norm_integ; % an opposite measurement
%
%   The number of columns of C, r, represents the number of vectors
%   from 0 (rays) in a n-dimensional space. [n,r]=size(C).
%   This averaged residual is calculated by averagin the squared 2-norm
%   distance from a point v in the dims-dimensional space to the convex
%   hull of the cone defined by the column vectors of C over the [0,1]^n
%   hypercube.
%
%   See also RANK, LSQNONNEG, INT_RES.

%   This file is part of AVERPOIN
%   (AVErage of squared Residuals for Positive-Input Networks)
%   $Revision: 1.5 $  $Date: 2016/25/01 20:19:00 $
%     Copyright (C) 2016  Richard R. Carrillo (University of Granada)
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
function res=int_res_num(cone_rays,n_pts_dim,plot_flag)
n_dims=size(cone_rays,1); % Num. of ray coordinates

if n_pts_dim^n_dims >= 1e6
   warning('This function is intended for testing and illustration of cones in a small dimensional space and achieves a limited precision. Reduce the number of rows of cone_rays input matrix or n_pts_dim to obtain the result in a reasonable time.')
end

res=total_int_num(cone_rays,n_pts_dim); % sum the integral value of the squared distance over each region simplex

if nargin>2 && plot_flag
   plot_cones_num(cone_rays,n_pts_dim) % Plot the regions in the hypercube
end
