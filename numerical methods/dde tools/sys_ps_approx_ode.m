function sys = sys_ps_approx_ode(sysd,tau,d,M)
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
%   tau: vector of time delays
%   d: system dimension
%   M: Chebyshev mesh dimension
% Output:
%   sys: linear ODE data structure
%    -> f: vector of drift terms @(t,p)
%    -> g: vector of diffusion terms @(t,p)

% Assure compatibility
if ~iscell(sysd.Fd)
    sysd.Fd = {sysd.Fd};
end
if ~iscell(sysd.Gd)
    sysd.Gd = {sysd.Gd};
end

% Spectral differentiation matrix
tau_max = max(tau);
Dm = 1/tau_max*cheb_diff(M);
Dm(end,:) = 0;
D0 = kron(Dm,eye(d));

% linear ODE system definition
sys.f = @(x,t,p) (D0 + psa_del_interp(sysd.F,sysd.Fd,tau,M,t,p))*x +...
    [zeros((M-1)*d,1); sysd.f(t,p)];
sys.g = @(x,t,p) psa_del_interp(sysd.G,sysd.Gd,tau,M,t,p)*x +...
    [zeros((M-1)*d,1); sysd.g(t,p)];

end
