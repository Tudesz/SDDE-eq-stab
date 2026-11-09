function [Q, R, iu, il] = vecmat_symm(n)
%VECMAT_SYMM Projection matricies to convert vectorized representations of
% symmetric matrices to the minimal, non-redundant form
% Input:
%   n: matrix dimension
% Output:
%   Q: projection matrix to unique coordinates R^{n^2} -> R^{(n^2+n)/2} 
%   R: projection matrix to redundant coordinates R^{(n^2+n)/2} -> R^{n^2} 
%   iu: indicies of upper diagonal elements
%   il: indicies of lower diagonal elements

% enumeration of matrix elements
Mi = reshape(1:n^2,n,n);
Mir = tril(Mi) + tril(Mi).' - diag(diag(Mi));
ij = reshape(Mir,[],1);
ir = unique(ij);
nr = length(ir);

% Transformation matricies
Q = zeros(nr,n^2);
for i = 1:nr
    Q(i,ir(i)) = 1;
end
R = zeros(n^2,nr);
for i = 1:n^2
    R(i,ij(i)==ir) = 1;
end

% Redundancy indicies
iu = reshape(triu(Mi)-diag(diag(Mi)),[],1);
iu(iu==0) = [];
il = reshape(tril(Mi)-diag(diag(Mi)),[],1);
il(il==0) = [];

end