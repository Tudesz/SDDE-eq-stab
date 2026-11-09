function [res,stab] = condfunc_cov_sim(X,Y,data)
%CONDFUNC_COV_SIM Second moment stability calculation of stochastic DDE
%equilibrium points via EM simulations
% Input: 
%   X: first varied system parameter
%   Y: second varied system parameter
%   data: data structure describing the stochastic differential equation
%    -> sys: data sturcture of the SODE:
%       -> f: drift term @(x,xd,t,p)
%       -> g: cell array of diffusion terms @(x,xd,t,p)
%    -> sysL: data sturcture of the linear SDDE:
%       -> F: matrix of drift terms @(t,p)
%       -> Fd: matrices of delayed drift terms {@(t,p)}
%       -> G: matrix of diffusion terms @(t,p)
%       -> Gd: matrices of delayed diffusion terms {@(t,p)}
%    -> p0: default parameter vector
%    -> pi_X: index of X in system parameter vector
%    -> pi_Y: index of Y in system parameter vector
%    -> sim: simulation options data structure for "sdde_EM_sim.m"
%    -> Xhist: initial history function
%    -> std_tol: std tolerance for allowing stability (default 2)
% Output:
%   res: output structure compatible with zd_trianglestab.m
%    -> aceig: ratio of the expected and simulated standard deviations
%    -> avg0: expected solution mean
%    -> std0: expected standard deviation
%    -> avg: simulated solution mean
%    -> std: simulated standard deviation
%   stab: stability flag (1 if res.aceig > data.std_tol)

% Default options
if ~isfield(data,'std_tol')
    data.std_tol = 2;
end

% Unpack data
par = data.p0;
par(data.pi_X) = X;
par(data.pi_Y) = Y;
T = data.sim.T;
F = data.sysL.F;
if ~iscell(data.sysL.Fd)
    Fd = {data.sysL.Fd};
else
    Fd = data.sysL.Fd;
end

% Expected stationary standard deviation
J0 = F(T,par); d = size(J0,1);
for i = 1:length(data.tau)
    J0 = J0 + Fd{i}(T,par);
end
x0 = -J0\data.sysL.f(T,par); % equilibrium point
x0t = @(t) repmat(x0,data.M,1);
JC = pde_cov_jac(T,data.sysL,par,data.tau,data.M);
bC = pde_cov_const(T,x0t,data.sysL,par,data.tau,data.M);
C0 = reshape(-JC\reshape(bC,[],1),data.M*d,data.M*d); 
c0 = C0(end-d+1:end,end-d+1:end); % second moment
s0 = sqrt(diag(c0-x0*x0.')); % standard deviations

% Simulate sytstem behaviour
[~,U,~,~] = sdde_EM_sim(data.sys,par,data.Xhist,data.tau,data.sim);
% figure(); plot(squeeze(U(1,:,:)))
Eu = zeros(d,1); Su = zeros(d,1);
for i = 1:d
    pUi = fitdist(squeeze(U(i,end,:)),'normal');
    Eu(i) = pUi.mu; Su(i) = pUi.sigma;
end

% Evaluate stability
if any(s0 == 0) || any(isnan(s0)) || any(imag(s0) ~=0)
    res.aceig = max(Su);
else
    res.aceig = max(Su./s0);
end
res.aceig = res.aceig - data.std_tol;
if res.aceig > 0 || norm(data.Xhist(0)) < norm(Eu)
    stab = 1; 
else
    stab = 0; % first and second moment stable point
end

% Output data
res.avg0 = x0;
res.std0 = s0;
res.avg = Eu;
res.std = Su;

end