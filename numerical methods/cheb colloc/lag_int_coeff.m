function L = lag_int_coeff(tm, tq)
%LAG_INT_COEFF Coefficients of Lagrange interpolation polynamials on mesh 
% tm, for querry points tq: f(tq) = L*f(tm)
% Input:
%   tm: default mesh points
%   tq: query point
% Output:
%   L: Lagrange interpolation coefficient matrix (p x q)

q = numel(tm);
p = numel(tq);

zi = repmat(reshape(tq, [p 1 1]), [1 q q]);
sj = repmat(reshape(tm, [1 q 1]), [p 1 q]);
sk = repmat(reshape(tm, [1 1 q]), [p q 1]);

t1 = zi-sk;
t2 = sj-sk;
idx = find(abs(t2)<=eps);
t1(idx) = 1;
t2(idx) = 1;

L = prod(t1./t2, 3);

end
