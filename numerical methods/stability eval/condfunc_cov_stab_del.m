function [res,stab] = condfunc_cov_stab_del(X,Y,data)
%CONDFUNC_COV_STAB Second moment stability calculation of stochastic DDE
%equilibrium points
% Input: 
%   X: first varied system parameter
%   Y: second varied system parameter
%   data: data structure describing the stochastic differential equation
%    -> sys: data sturcture of the linear SDDE:
%       -> F: matrix of drift terms @(t,p)
%       -> Fd: matrices of delayed drift terms {@(t,p)}
%       -> G: matrix of diffusion terms @(t,p)
%       -> Gd: matrices of delayed diffusion terms {@(t,p)}
%    -> tau: vector of time delays (or a function in parameters)
%    -> M: dimension of the Chebyshev mesh used for discretization
%    -> p0: default parameter vector
%    -> pi_X: index of X in system parameter vector
%    -> pi_Y: index of Y in system parameter vector
%    -> Q: left projection matrix to symmetrical part (optional)
%    -> R: right projection matrix to symmetrical part (optional)
% Output:
%   res: output structure compatible with zd_trianglestab.m
%    -> aceig: maximal real part of charachteristic eigenvalues
%    -> lam_c: critical eigenvalue
%   stab: stability flag (1 if res.aceig > 0)

% Unpack data
par = data.p0;
par(data.pi_X) = X;
par(data.pi_Y) = Y;
if isnumeric(data.tau)
    tau = data.tau;
else
    tau = data.tau(par);
end

% Evaluate stability
JC = pde_cov_jac(0,data.sys,par,tau,data.M);
if isfield(data,'Q')
    lam = eig(data.Q*JC*data.R);
else
    lam = eig(JC);
end
[~,ic] = max(real(lam));
res.aceig = real(lam(ic));
res.lam_c = lam(ic);
if res.aceig > 0
    stab = 1;
else
    stab = 0;
end

end