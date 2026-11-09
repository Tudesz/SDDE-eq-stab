function eq = sdde_char_sn(sys,tau,p0)
%SDDE_CHAR_EQ Evaluate the characteristic equation of the correlation
% boundary value problem of a linear stochastic dde considering the 
% limiting case of lambda = 0
% TODO: account for multiple delay terms and noise sources
% TODO: account for explicit time dependence in the functions found in sys
% Input: 
%   sys: data sturcture of the linear SDDE:
%       -> F: matrix of drift terms @(t,p)
%       -> Fd: matrices of delayed drift terms {@(t,p)}
%       -> G: matrix of diffusion terms @(t,p)
%       -> Gd: matrices of delayed diffusion terms {@(t,p)}
%   tau: vector of time delays
%   p0: parameter vector
% Output:
%   eq: evaluation of the characteristic equation det(M(lambda))=0

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

% Construct the system matricies
Ai = tau/2*[-ah -bhP; bhP ah]; %Int(A(lambda,delta))
b =[aht+be2 bht]; % b^T(lambda)

% Account for singularities considering the symmetry of f(tau/2)
[Q, ~, iu, il] = vecmat_symm(d); % projection to the minimum coordinates necessary
I_f = eye(d^2);
bns = [Q*b;  I_f(iu,:)-I_f(il,:) zeros(length(iu),d^2)]; % singular conditions replaced

% Characteristic equation
M = [bns*expm(Ai); eye(d^2) -eye(d^2)];
eq = det(M);

end