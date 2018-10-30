%CONVHULLN0 N-D Convex hull of a polytope containing the origin vertex.
%   [K,V] = CONVHULLN0(X) returns the indices K of the points in X that
%   comprise the facets of the convex hull of X.
%   X is an n-by-m array representing m points in n-D space. The origin
%   point (0) is not explicitly specified but it is assumed to be included.
%   If the convex hull has p facets then K is d-by-p. This matrix contains
%   indices to the points (columns) in x. The index to the first vertex
%   in X is 1. The index 0 referres to the origin vertex.
%   The main difference with respect convhulln is that this function works
%   even when num. of non-collinear points is less than the num. of
%   dimensions (by performing a previous dimension reduction).
%   If this is the case, d could be lower than n, otherwise, d=n.
%   V returns the polytope volume.
function [k0,v] = convhulln0(x,options)
space_ndims=size(x,1); % num. of space dimensions
indep_rows=rank(x); % number of linearly independent rows (rays)
% convhulln works if the number of linearly independent points is higher
% than the space dimension
if indep_rows < space_ndims % we have less than n+1 collinear points: we reduce dimensionality
   if indep_rows ~= 0
      % Calculate an orthonormal basis for the space which the cone spans
      orth_basis=orth(x); % SVD performed for rank() could have been used for orth() for effciency purposes
      x_red=(x'*orth_basis)'; % Update coordinates of points
   else
      x_red=[]; % No dimension reduction has been performed
   end
   dim_red=true; % Note down that we have performe dim. reduction
else
   x_red=x; % No dimension reduction has been performed
   dim_red=false;
end

% cone_ndims is the number of space dimensions that the cone actually
% spans. If the input cone rays define a full-dimenisional cone,
% cone_ndims is equal to space_ndims. cone_ndims is a natural number of
% the interval [0, space_ndims].
cone_ndims=size(x_red,1);

% In dimension n conhulln requires n linearly-independent
% points plus another point which is different from the others
% (reference), in our case the reference is the origin
% (x_red_0 is the total point set).
% We take the points defined by these vectors and append the origin
% point. The resultant matrix has enough non-collinar points to
% calculate the convex hull of the corresponding polytope (through
% the convhull function).
x_red_0=[zeros(cone_ndims,1) x_red];

% convhulln does not work for one-dimensional space, so check it
if cone_ndims>1
   % In this software it is the convention to store the vectors in the
   % columns of the matrices. However, convhull expect that the input
   % points are encoded by the input matrix rows.
   if nargin > 1
      opt_par=options;
   else
      opt_par=[]; % default options (see convhulln())
   end
   [k_trans,v]=convhulln(x_red_0',opt_par);
   k=k_trans'; % for convhulln each row of k is a facet. For conhulln0 each column is a facet
else % Compute convex hull for one-dimensional space
   if size(x_red_0,2)>1
      [max_coord,max_ind]=max(x_red_0);
      [min_coord,min_ind]=min(x_red_0);
      k=[min_ind max_ind];
      v=max_coord-min_coord;
   else % If all the point have the same coordinates, x becomes an empty matrix
      k=[1];
      v=0;
   end
end
k0=k-1; % Make 0 the index to the origin point
if dim_red % The input polytope is not full dimensional
   v=0;
end