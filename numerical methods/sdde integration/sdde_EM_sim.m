function [t,U,Ud,dWt] = sdde_EM_sim(sys,p0,X0,tau,sim)
%SDDE_EM_SIM Euler-Maruyama simulation routine with fixed Wiener
% process for stochastic DDEs given in the form: 
% dX(t) = f(X(t),X(t-tau),p,t)dt + g(X(t),X(T-tau),p,t)dW(t)
% Input:
%   sys: data sturcture of the SODE:
%    -> f: drift term @(x,xd,t,p)
%    -> {g}: cell array of diffusion terms @(x,xd,t,p)
%   p0: system parameter vector
%   X0: initial condition / history function @(t) x [-tau_max 0] 
%      (overwritten if sim.X0 is available)
%   tau: fixed, point time delays
%   sim: simulation option data structure
%    -> T: end time of the simulation
%    -> R: Brownian samples within an intagration step (default 1)
%    -> dt: solver timestep (default 1e-3, overwritten if dW is available)
%    -> n: number of trajectories to simulate (default 1, 
%       overwritten if dW is available)
%    -> X0: cell array of pre-generated initial functions (optional)
%    -> logs: log integration progress (default false)
% Output:
%   t: time signal of the solution
%   U: ensamble of generated stochastic solution signals
%   Ud: ensamble of generated stochastic delayed term signals
%   dWt: realisations of the Wiener increments on t

% Default settings, options
if ~isfield(sim,'dt')
    sim.dt = 1e-3;
end
if ~isfield(sim,'n')
    sim.n = 1;
end
if ~isfield(sim,'R')
    sim.R = 1;
end
if ~isfield(sim,'logs')
    sim.logs = false;
end
if ~iscell(sys.g)
    sys.g = {sys.g};
end

if isfield(sim,'X0')
    X0 = @(t) sim_hist(t,sim.X0);
end

% Noise process and initial condition initialization
[t, ~, dW] = brownian_path(sim.T,sim.dt,sim.n,length(sys.g));

% Forward time integration
[t, U, Ud, dWt] = sdde_em(sys.f,sys.g,p0,X0,tau,t,dW,sim);

end

% Auxiliary functions
function xt = sim_hist(t,X0)
%SIM_HIST unpack pre-generated hsitory function cell array
% Input:
%   t: querry time (n)
%   X: cell array of history functions {n}
% Output:
%   xt: evaluations of the state at t [d x n]

d = size(X0{1}(0),1);
n = length(t);
xt = zeros(d,n);
for i = 1:n
    xt(:,i) = X0{i}(t(i));
end

end
