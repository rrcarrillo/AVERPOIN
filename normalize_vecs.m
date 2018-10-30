% col_nvecs=normalize_vecs(col_vecs) normalizes the column vectors of
% col_vecs.
% This function is equivalent to col_nvecs=normc(col_vecs) except that
% normalize_vecs returns a zero vector for each column vector which is 0.
function col_nvecs=normalize_vecs(col_vecs)
% Calculate the norm of column vectors
vec_norms=sqrt(sum(col_vecs.^2,1));
% Avoid norms of length 0
vec_norms(vec_norms==0)=1;
% Divide vectors by their normals
col_nvecs=bsxfun(@rdivide,col_vecs,vec_norms);
