function [th, phi] = corr_stat_asol(sys,tau,p0,compact)
%CORR_STAT_ASOL Return the analytic solution of the covariance function
%for scalar stochastic DDEs (only solving the homogenous problem!)
% TODO: account for multiple delay terms and noise sources
% TODO: account for explicit time dependence in the functions found in sys
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
%   compact: if true only return the covariance at th=0 (default false)
% Output:
%   th: time mesh to evaluate phi
%   phi: analytic solution of E[x(t)x(t+theta)] 

if nargin < 4
    compact = false;
end
if ~iscell(sys.Fd)
    sys.Fd = {sys.Fd};
end
if ~iscell(sys.Gd)
    sys.Gd = {sys.Gd};
end

% Unpack function definitions
F = sys.F(0,p0);
Fd = sys.Fd{1}(0,p0);
ft = sys.f(0,p0); 
G = sys.G(0,p0);
Gd = sys.Gd{1}(0,p0);
gt = sys.g(0,p0);
tau = tau(1); % only one delay is considered
d = size(F,1);

% Diagonalize the governing equations: phi'(th) = -ah*phi(th) - bh*phi(-tau-th)
P = vecmat_transp(d); % permutation matrix for transposing
ah = kron(eye(d),F);
bh = kron(eye(d),Fd)*P;
M = [-ah -bh; bh ah]; % fg' = M*fg
[V,D] = eig(M); % fg(t) = A*expm(M*t)
lam = diag(D); % characteristic exponents
nl = length(lam)/2;

% Stationary solution
x0 = -(F+Fd)\ft; % (a + b)*x0 = -c
gh = kron(eye(d),gt)*gt + (kron(eye(d),ft) + kron(ft,eye(d)) + ...
    kron(eye(d),G)*kron(gt,eye(d)) + kron(G,eye(d))*kron(eye(d),gt) + ...
    kron(eye(d),Gd)*kron(gt,eye(d)) + kron(Gd,eye(d))*kron(eye(d),gt))*x0;
fg0 = M\kron([1;-1],kron(eye(d),ft)*x0); % [cx -cx].' = M*fg0
f0 = fg0(1:nl);
g0 = fg0(nl+1:end);

% Solve for the coefficents based on the boundary conditions
vF = V(1:nl,:);
vG = V(nl+1:end,:);
fb = (kron(eye(d),F) + kron(F,eye(d)) + kron(eye(d),G)*kron(G,eye(d)) + ...
    kron(eye(d),Gd)*kron(Gd,eye(d))); % 2a + alpha^2 + beta^2
gb = (kron(eye(d),Fd)*P + kron(Fd,eye(d)) + kron(eye(d),G)*kron(Gd,eye(d)) ...
    + kron(eye(d),Gd)*kron(G,eye(d))*P); % 2b + 2*alpha*beta
Bc1 = vF-vG; % f(0) = g(0)
Bc2 = exp(lam.'*tau/2).*(fb*vF+gb*vG); % fb*f(tau/2) + gb*g(tau)/2 
bc = -gh - fb*f0 - gb*g0; % eta^2 -fb*f0 - gb*g0

% Account for singularities
[Q, ~, iu, il] = vecmat_symm(d); % projection to the minimum coordinates necessary
Bc2 = [Q*Bc2; exp(lam.'*tau/2).*(vF(iu,:)-vF(il,:))]; % phi(0) = phi(0).'
c = [Bc1; Bc2]\[zeros(nl,1); Q*bc; zeros(length(iu),1)];

% Reconstruct solution
if ~compact
    f = @(del) f0 + real(c.'.*exp(lam.'*del)*vF.').';
    g = @(del) g0 + real(c.'.*exp(lam.'*del)*vG.').';
    delta = linspace(0,tau/2,1000);
    fd = zeros(d^2,length(delta));
    gd = zeros(d^2,length(delta));
    for i = 1:length(delta)
        fd(:,i) = f(delta(i)); 
        gd(:,i) = g(delta(i)); 
    end
    % figure(); plot(delta,fd,-delta,gd)
    th = [-flip(delta) delta(2:end)] - tau/2;
    phi = [flip(gd,2) fd(:,2:end)];
else
    th = 0; 
    phi = reshape(f0 + real(c.'.*exp(lam.'*tau/2)*vF.').',d,d);
end

% Check boundary conditions
% disp(norm(f(0)-g(0)))
% disp(norm(fb*f(tau/2)+gb*g(tau/2)-bc))
% df = @(del) (lam.'.*c.'.*exp(lam.'*del)*vF.').';
% dg = @(del) (lam.'.*c.'.*exp(lam.'*del)*vG.').';
% disp(norm([df(tau/2);dg(tau/2)]-M*[f(tau/2);g(tau/2)]))

end