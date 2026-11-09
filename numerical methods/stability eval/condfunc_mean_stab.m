function [res,stab] = condfunc_mean_stab(X,Y,data)
%CONDFUNC_MEAN_STAB First moment stability calculation of stochastic ODE
%equilibrium points
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
%    -> aceig: maximal real part of charachteristic eigenvalues
%    -> lam_c: critical eigenvalue
%   stab: stability flag (1 if res.aceig > 0)

% Unpack data
par = data.p0;
par(data.pi_X) = X;
par(data.pi_Y) = Y;

% Evaluate stability
A = data.sys.F(0,par);
lam = eig(A);
[~,ic] = max(real(lam));
res.aceig = real(lam(ic));
res.lam_c = lam(ic);
if res.aceig > 0
    stab = 1;
else
    stab = 0;
end

end