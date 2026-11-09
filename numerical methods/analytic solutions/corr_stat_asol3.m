function [A, f0] = corr_stat_asol3(sys,tau,p0)
%CORR_STAT_ASOL2 Return the analytic solution of the covariance function
%for scalar stochastic DDEs (using matrix exponential formalism)
% Input:
%   sys: data sturcture of the linear SDDE:
%    -> F: matrix of drift terms @(t,p)
%    -> Fd: matrices of delayed drift terms {@(t,p)}
%    -> f: vector of inhomogenous drift terms @(t,p)
%    -> G: matrix of diffusion terms @(t,p)
%    -> Gd: matrices of delayed diffusion terms {@(t,p)}
%    -> g: vector of inhomogenous diffusion terms @(t,p)
%   tau: vector of time delays
%   p0: system parameter vector
% Output:
%   A: coeffients in the matrix exponential
%   f0: correlation function at 0

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
gt = sys.g(0,p0);
tau = tau(1); % only one delay is considered
d = size(F,1);

% Matrix form of the governing BVP
P = vecmat_transp(d); % permutation matrix for transposing
[Q, ~, iu, ~] = vecmat_symm(d); % projection to the minimum coordinates necessary
Qc = eye(d^2); Qc = Qc(iu,:); % complement of the projection matrix
ah = kron(eye(d),F); % I otimes a
bh = kron(eye(d),Fd)*P; % I otimes b and transpose
A = [-ah -bh; bh ah]; % fg' = A*fg
B = [kron(eye(d),F) + kron(F,eye(d)) + kron(G,G) + kron(Gd,Gd) ...
    kron(eye(d),Fd)*P + kron(Fd,eye(d)) + kron(Gd,G) + kron(G,Gd)*P];

% Stationary solution
M = [Q*B; Qc*[P-eye(d^2) zeros(d^2)]]*expm(A*tau/2)*[eye(d^2); eye(d^2)];
f0 = M\[-Q*reshape(gt*gt.',[],1); zeros(d*(d-1)/2,1)];

end