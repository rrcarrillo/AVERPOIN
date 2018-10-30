% This function is used to determine if each vector of a set is pointing 
% towards the outer side of a cone facet
% outer_f = facet_orient(cone_center, facet_u, vec)
% cone_center is a column vector with the coordinates of an approximate
% cone center ray (a ray that is inside the cone)
% facet_u are the vertices (rays) defining the facet of the cone to
% consider.
% pts is a matrix. Each column represents a point whose position with regard
% to the facet orientation is calculated.
% outer_f is a vector. Each of its elements corresponds to a point (column
% of pts). If outer_f element is true, the point is facing the outer side of the
% facet facet_u . If it is false, the point is facing the inner side of the facet.
function outer_f=facet_orient(cone_center,facet_u,vecs)
large_eps=eps(1000); % Values smaller than this are considered zero

% Normalize center ray (in order to calculate correctly the distances to the points)
cone_center_norm=normalize_vecs(cone_center);

% Project center ray on the current facet
base_facet_u=gschmidt(facet_u); % Calculate the cone facet subspace othonormal basis
% cone_center_on_facet is a column vector containing the coordinates of the
% ray resulting from projecting cone_center ray on the facet defined by facet_u.
cone_center_on_facet=normalize_vecs(((cone_center_norm'*base_facet_u)*base_facet_u')');

% Function to compare two vectors using absolute value tolerance
isequal_tol = @(x,y) all(abs(x-y) < large_eps);

% If the position of pcenter_u do not change after proyecting it on the
% facet, the cone has no volume (it is not full dimentional)
full_dim_cone=~isequal_tol(cone_center_norm, cone_center_on_facet);

% If the cone has no volume (it is not full dimentional), use the
% facet vertex order the determine facet orientation to the points in pts
% Put each column in a cell
pts_cell=num2cell(vecs,1);
% For each cell of pts_cell (a point) append it to the facet_u matrix as
% a new column and calculate the determinat of the resultant matrix
pts_orient_det=cellfun(@(pt) det([facet_u pt]),pts_cell,'UniformOutput',true);
% Check the sign of the calcualted determinants to know the orientation
% with respect the point
outer_f=pts_orient_det<0;

if full_dim_cone % If the cone has volume we can ensure that the facer vertices are properly ordered
   if det([facet_u cone_center_norm]) < 0
      % Ups! The facet vertices are not properly ordered to define the
      % facet orientation (as convhull sometimes does)
      outer_f=~outer_f; % Negate the function output
   end
else
   warning('Non-full-dimensional: facet normal direction may be incorrectly calculated')
end
