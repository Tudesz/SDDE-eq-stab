function A = dde_fd_igener(x,xd,p,M,sys)
%DDE_FP_IGENER Infinitesimal generator of a locally linearized neutral
% delay differential equation (modified version of ndde_inf_gener using 
% a finite difference scheme)
% Input:
%   x: solution to linearize about (current state)
%   xd: solution to linearize about (delayed states)
%   p: parameter vector
%   M: Chebyshev mesh dimension
%   sys: system definition data structure
%    -> f: vector field and its Jacobians
%    -> tau: time delay and its parameter Jacobian
%    -> tau_no: number of distinct time delays
% Output:
%   A: infinitesimal generator of the local linear system

% unpack time delays
nt = sys.tau_no;
tau = zeros(1,nt);
for i = 1:nt
    tau(i) = feval(sys.tau,p,i,2);
end
tau_max = max(tau);
tM = linspace(-tau_max,0,M);

% initialize the A matrix
n = size(x,1);
Dm = (M-1)/tau_max*(diag(ones(M-1,1),1)-eye(M));
A = kron(Dm,eye(n));

% Function definition for vector space Jacobians
Ji = @(i) feval(sys.f,x,xd,p,2,i); 

% Fill up the A matrix
A(end-n+1:end,:) = 0;               % clear rows corresponding to the present
A(end-n+1:end,end-n+1:end) = Ji(0); % effect of non delayed terms
for i = 1:nt
    lq = lag_int_coeff(tM,-tau(i));
    A(end-n+1:end,:) = A(end-n+1:end,:) + kron(lq,Ji(i));
end


end