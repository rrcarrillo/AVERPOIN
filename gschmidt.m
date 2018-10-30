% [Q,R]=gschmidt(A)
% Gram-schmidt orthogonalization
% Input: A is an m by n matrix of full rank m<=n whose columns are vectors
% Output: an m-by-n upper triangular matrix (R) with the inner products
% q_i'*a_j
% and an m-by-m unitary matrix of orthonormal vectors (the columns of Q)
% so that A = Q*R.
function [Q,R]=gschmidt(A)
[m n] = size(A);
Q = zeros(m,n);
R = zeros(n);
for j=1:n
    v=A(:,j); % v begins as column j of A
    for i=1:j-1
        R(i,j)=Q(:,i)'*A(:,j); % modify A(:,j) to v for more accuracy
        v=v-R(i,j)*Q(:,i); % substract the proyection
    end % v is now perpendicular to all of q_1,...,q_j-1
    R(j,j)= norm(v);
    if R(j,j)~=0
        Q(:,j)=v/R(j,j); % normalize v (if possible) to be the next unit vector q_j
    else
        Q(:,j)=v;
    end
end