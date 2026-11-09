function P = vecmat_transp(n)
%VECMAT_TRANSP Calulate the permutation matrix for transposing a matrix in
%vectorized form vec(A.') = P*vec(A)
% Input:
%   n: matrix dimension
% Output:
%   P: permutation matrix (n^2 by n^2)

ij = reshape(reshape(1:n^2,n,n).',[],1);
P = zeros(n^2);
for i = 1:n^2
    P(i,ij(i)) = 1;
end

end