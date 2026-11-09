function sys = sys_ps_approx_lode(sysd,tau,d,M)
%SYS_PS_APPROX_LODE Pseudo-spactral ODE approximation of a linear DDE with
% fixed delays
% Input:
%   sysd: linear delay differential equation data structure
%    -> F: matrix of drift terms @(t,p)
%    -> Fd: matrices of delayed drift terms {@(t,p)}
%    -> f: vector of inhomogenous drift terms @(t,p)
%    -> G: matrix of diffusion terms @(t,p)
%    -> Gd: matrices of delayed diffusion terms {@(t,p)}
%    -> g: vector of inhomogenous diffusion terms @(t,p)
%   tau: vector of time delays (or a function of parameters)
%   d: system dimension
%   M: Chebyshev mesh dimension
% Output:
%   sys: linear ODE data structure
%    -> F: matrix of drift terms @(t,p)
%    -> f: vector of inhomogenous drift terms @(t,p)
%    -> G: matrix of diffusion terms @(t,p)
%    -> g: vector of inhomogenous diffusion terms @(t,p)

% Assure compatibility
if ~iscell(sysd.Fd)
    sysd.Fd = {sysd.Fd};
end
if ~iscell(sysd.Gd)
    sysd.Gd = {sysd.Gd};
end
if isnumeric(tau)
    tau = @(p) tau;
end

% Spectral differentiation matrix
Dm = cheb_diff(M);
Dm(end,:) = 0;
D0 = kron(Dm,eye(d));

% linear ODE system definition
sys.F = @(t,p) 1/max(tau(p))*D0 + ...
    psa_del_interp(sysd.F,sysd.Fd,tau(p),M,t,p);
sys.f = @(t,p) [zeros((M-1)*d,1); sysd.f(t,p)];
sys.G = @(t,p) psa_del_interp(sysd.G,sysd.Gd,tau(p),M,t,p);
sys.g = @(t,p) [zeros((M-1)*d,1); sysd.g(t,p)];

end
