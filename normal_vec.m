% normal_vec(u) return a vector normal to the hyperplane defined by the
% column vectors of u (rays).
% The normal (perpendicular) vector is calculated by means of the Laplace
% expansion of the determinant:
% https://ef.gy/linear-algebra:normal-vectors-in-higher-dimensional-spaces
function n=normal_vec(u)
[n_dims,n_vectors]=size(u);
if n_dims~=n_vectors+1
   error('The number of columns of u must be equal to its number of rows minus 1')
end
n=zeros(n_dims,1);
for n_dim=1:n_dims
   sign_coeff=(-1)^(n_dims+n_dim);
   cofactor=det(u((1:n_dims)~=n_dim,:))*sign_coeff;
   n(n_dim)=cofactor;
end
% If n_dims is large, the length of the resultant vector may be very small,
% so normalized n to a length of one
n=normalize_vecs(n);
