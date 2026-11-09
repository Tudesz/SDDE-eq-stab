function dX = sdde_ps_approx(X,t,p,f,d,M,tau)
%SDDE_PS_APPROX Discretization of a stochastic delay differential equations 
% exploiting a pseudo-spectral approximation 
% Input:
%   X: discretized solution on the whole solution tail (d*M x d)
%   t: current time
%   p: parameter vector
%   d: state dimension
%   f: function handle for the deterministic part of the DDE rhs (x,xd,t,p)
%   M: Chebyshev mesh dimesnion used for discretizationS
%   tau: fix, point time delays present in the system
% Output:
%   dX: approximate derivatives for the whole discretized solution tail

% unpack time delays
nt = length(tau);
tau_max = max(tau);

% Spectral differentiation matrix
Dm = 1/tau_max*cheb_diff(M);
D = kron(Dm,eye(d));
D(end-d+1:end,:) = [];

% Evaluate the delayed terms
n = size(X,2);
x = X(end-d+1:end,:);
xd = zeros(d,n,nt);
for i = 1:nt
    tq = 1-tau(i)/tau_max; % rescaled to [0 1]
    lq = lag_coeff(M,tq);
    xd(:,:,i) = kron(lq,eye(d))*X;
end

% Evaluate the extended ODE rhs
dX = [D*X; f(x,xd,t,p)];

end