function dX = sdde_fd_approx(X,t,p,f,d,M,tau)
%SDDE_FD_APPROX Discretization of a stochastic delay differential equations 
% exploiting a finite difference method 
% Input:
%   X: discretized solution on the whole solution tail (d*M x n)
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
tM = linspace(-tau_max,0,M);

% Forward differentiation matrix
Dm = (M-1)/tau_max*(diag(ones(M-1,1),1)-eye(M));
D = kron(Dm,eye(d));
D(end-d+1:end,:) = [];

% Evaluate the delayed terms
n = size(X,2);
x = X(end-d+1:end,:);
xd = zeros(d,n,nt);
for i = 1:nt
    lq = lag_int_coeff(tM,-tau(i));
    xd(:,:,i) = kron(lq,eye(d))*X;
end

% Evaluate the extended ODE rhs
dX = [D*X; f(x,xd,t,p)];

end