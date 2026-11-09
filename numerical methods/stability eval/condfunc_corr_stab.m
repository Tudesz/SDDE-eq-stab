function [res,stab] = condfunc_corr_stab(X,Y,data)
%CONDFUNC_CORR_STAB Second moment stability calculation of stochastic DDE
%equilibrium points based on positive definiteness of stationary
%correlation functions
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
Jphi = pde_corr_jac(0,data.sys,par,tau,data.M);
res.aceig = det(Jphi);
if res.aceig > 0
    stab = 1;
else
    stab = 0;
end

end