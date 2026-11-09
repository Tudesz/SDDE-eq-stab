function F = psa_del_interp(A,B,tau,M,t,p)
%PSA_DEL_INTERP: Form the pseudo-spectral approximation of the DDE rhs
% Input:
%   A: coefficient matrix of x(t)
%   B: coefficient matrices of x(t-tau_i)
%   tau: time delay vector
%   M: Chebyshev mesh resolution
% Output:
%   F: M*d by M*d coefficient matrix of the discretized history

% Output initialization
d = size(A(t,p),1);
F = zeros(d*M);

% Fill up the last rows
F(end-d+1:end,end-d+1:end) = A(t,p);
tau_max = max(tau);
for i = 1:length(tau)
    tq = 1-tau(i)/tau_max; % rescaled to [0 1]
    lq = lag_coeff(M,tq); % Lagrange coefficients
    F(end-d+1:end,:) = F(end-d+1:end,:) + kron(lq,B{i}(t,p));
end

end