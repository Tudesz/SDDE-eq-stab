function phi = corr_scalar_asol(sys,tau,p0)
%CORR_SCALAR_ASOL Return the analytic solution of the covariance function
%for scalar stochastic DDEs
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
% Output:
%   phi: analytic solution of E[x(t)x(t+theta)] 

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

% Inhomogenous terms
x0 = - ft/(F+Fd);
eta = sqrt(gt^2 + (2*ft + 2*G*gt + 2*Gd*gt)*x0);

% Analytic solution (homogenous)
% phi = @(th)-(eta^2*(sqrt(F^2 - Fd^2)*cosh(sqrt(F^2 - Fd^2)*(th + tau/2)) ...
%     - (F + Fd)*sinh(sqrt(F^2 - Fd^2)*(th + tau/2))))/(sqrt(F^2 - Fd^2) ...
%     *(2*F + 2*Fd + (G + Gd)^2)*cosh(1/2*sqrt(F^2 - Fd^2)*tau) ...
%     - (F + Fd)*(2*F - 2*Fd + (G - Gd)^2)*sinh(1/2*sqrt(F^2 - Fd^2)*tau));

% Analytic solution (non-homogenous)
phi = @(th) (F.^2+(-1).*Fd.^2).^(-1/2).*((F.^2+(-1).*Fd.^2).^(1/2).*(2.*F+2.* ...
  Fd+(G+Gd).^2).*cosh((1/2).*(F.^2+(-1).*Fd.^2).^(1/2).*tau)+(-1).*( ...
  F+Fd).*(2.*F+(-2).*Fd+(G+(-1).*Gd).^2).*sinh((1/2).*(F.^2+(-1).* ...
  Fd.^2).^(1/2).*tau)).^(-1).*((F+(-1).*Fd).*((-1).*eta.^2.*(F+Fd)+ ...
  ft.*(2.*F+2.*Fd+(G+Gd).^2).*x0).*cosh((F.^2+(-1).*Fd.^2).^(1/2).*( ...
  (1/2).*tau+th))+ft.*x0.*((-1).*(F+(-1).*Fd).*(2.*F+2.*Fd+(G+Gd) ...
  .^2).*cosh((1/2).*(F.^2+(-1).*Fd.^2).^(1/2).*tau)+(F.^2+(-1).* ...
  Fd.^2).^(1/2).*(2.*F+(-2).*Fd+(G+(-1).*Gd).^2).*sinh((1/2).*(F.^2+ ...
  (-1).*Fd.^2).^(1/2).*tau))+(-1).*(F.^2+(-1).*Fd.^2).^(1/2).*(ft.*( ...
  G+Gd).^2.*x0+F.*((-1).*eta.^2+2.*ft.*x0)+Fd.*((-1).*eta.^2+2.*ft.* ...
  x0)).*sinh((F.^2+(-1).*Fd.^2).^(1/2).*((1/2).*tau+th)));

end