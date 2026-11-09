function JC = ode_cov_jac(t,sys,p0)
%ODE_COV_JAC Evaluate the Jacobian of the second moment equation that
%describes the linear stochastic ODE
% dx(t) = (Ft*x(t) + f)*dt +  (Gt*x(t) + gt)*dWt
% Input:
%   t: integration time
%   sys: data sturcture of the linear SODE:
%    -> F: matrix of drift terms @(t,p)
%    -> G: matrix of diffusion terms @(t,p)
%   p0: system parameter vector
% Output:
%   JC: Jacobian of the second moment equation wrt C

% Unpack function definitions
Ft = sys.F(t,p0);
Gt = sys.G(t,p0);
d = size(Ft,1);

% Differenatiation of the moment matrix
JC = kron(eye(d),Ft) + kron(Ft,eye(d)) +... % F*C + C*F.'
    kron(eye(d),Gt)*kron(Gt,eye(d)); % G*C*G.'

end