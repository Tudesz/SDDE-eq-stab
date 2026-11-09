function [res,stab] = condfunc_stat_stab(X,Y,data)
%CONDFUNC_MEAN_STAB Second-moment stability calculation of stochastic ODE
%equilibrium points based on positive definiteness of stationary
%correlation functions
% Input: 
%   X: first varied system parameter
%   Y: second varied system parameter
%   data: data structure describing the stochastic differential equation
%    -> sys: linear ODE data structure
%       -> F: matrix of drift terms @(t,p)
%    -> p0: default parameter vector
%    -> pi_X: index of X in system parameter vector
%    -> pi_Y: index of Y in system parameter vector
% Output:
%   res: output structure compatible with zd_trianglestab.m
%    -> aceig: determinant of the coefficient matrix (Psi)
%    -> feig: critical eigenvalue of the stationary covariance
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

% Estimate second-moment stability
[Psi, f0] = sdde_char_mat(data.sys,tau,par);
res.aceig = det(Psi);
res.feig = min(real(eig(f0)));
if res.aceig > 0
    stab = 1;
else
    stab = 0;
end

end