function Jc = corr_pert_jac(lambda,sys,tau,p0,M)
%CORR_PERT_JAC Jacobian of the linear correlation pertrubation BVP for
%stochastic delay differential equations
% Input:
%   lambda: characteristic multiplier
%   sys: data sturcture of the linear SDDE:
%       -> F: matrix of drift terms @(t,p)
%       -> Fd: matrices of delayed drift terms {@(t,p)}
%       -> G: matrix of diffusion terms @(t,p)
%       -> Gd: matrices of delayed diffusion terms {@(t,p)}
%   tau: vector of time delays
%   p0: parameter vector
%   M: discretization dimension
% Output:
%   Jc: Jacobian of the correlation function boundary value problem

if ~iscell(sys.Fd)
    sys.Fd = {sys.Fd};
end
if ~iscell(sys.Gd)
    sys.Gd = {sys.Gd};
end

% Unpack function definitions
F = sys.F(0,p0);
Fd = sys.Fd{1}(0,p0);
G = sys.G(0,p0);
Gd = sys.Gd{1}(0,p0);
tau = tau(1); % only one delay is considered
d = size(F,1);

% Vectorize the governing equations
P = vecmat_transp(d); % permutation matrix for transposing
ah = kron(eye(d),F);
bhP = kron(eye(d),Fd)*P;
aht = kron(eye(d),F) + kron(F,eye(d)) + kron(eye(d),G)*kron(G,eye(d));  % 2a + alpha^2 
be2 = kron(eye(d),Gd)*kron(Gd,eye(d)); % beta^2
bht = (kron(eye(d),Fd)*P + kron(Fd,eye(d)) + kron(eye(d),G)*kron(Gd,eye(d)) ...
    + kron(eye(d),Gd)*kron(G,eye(d))*P); % 2b + 2*alpha*beta

% Discretization of the BVP
theta = cheb_mesh(M,[-tau,0]);
Dm = cheb_diff(M);
Dm(end,:) = 0;
D = 1/tau*kron(Dm,eye(d^2));

% Inner points
A = zeros(d^2*M);
for i = 1:M-1
    ii = (i-1)*d^2+1:i*d^2; % indicies of phi(theta)
    A(ii,ii) = A(ii,ii) + lambda*eye(d^2) - ah;
    tq_i = -theta(i)/tau;
    fi = lag_coeff(M,tq_i); 
    A(ii,:) = A(ii,:) - kron(fi,bhP*exp(lambda*theta(i)));
end

% Boundary condtion
i0 = (M-1)*d^2+1:M*d^2;
A(i0,i0) = aht + be2*exp(-lambda*tau) - lambda*eye(d^2);
A(i0,1:d^2) = bht;

% Avoid redundant equations
[Q, ~, iu, il] = vecmat_symm(d); nr = d^2-size(Q,1);
A(i0,:) = [Q*A(i0,:); zeros(nr,d^2*M)];
A(end-nr+1:end,i0(iu)) = 1;
A(end-nr+1:end,i0(il)) = -1;

% Combine the equations
Jc = D - A;

end