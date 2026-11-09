function A = dde_lin_igener(sys,tau,M)
%DDE_LIN_IGENER Infinitesimal generator of a linear DDE
% Input:
%   sys: data sturcture of the linear SDDE:
%    -> F: matrix of drift terms
%    -> Fd: matrices of delayed drift terms
%   tau: vector of time delays
%   M: Chebyshev mesh dimension
% Output:
%   A: infinitesimal generator of the local linear system

% unpack time delays
nt = length(tau);
tau_max = max(tau);
tM = linspace(-tau_max,0,M);

% initialize the A matrix
n = size(sys.F,1);
Dm = (M-1)/tau_max*(diag(ones(M-1,1),1)-eye(M));
A = kron(Dm,eye(n));

% Fill up the A matrix
A(end-n+1:end,:) = 0; % clear rows corresponding to the present
A(end-n+1:end,end-n+1:end) = sys.F; % effect of non delayed terms
for i = 1:nt
    lq = lag_int_coeff(tM,-tau(i));
    A(end-n+1:end,:) = A(end-n+1:end,:) + kron(lq,sys.Fd{i});
end


end