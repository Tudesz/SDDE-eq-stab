clear; close all; clc;
%% Simulate the covariance of a linear, scalar SDDE
% dx(t) = (a*x(t)+b*x(t-tau) + c)*dt + 
%   (alpha*x(t)+beta*x(t-tau) + gamma)*dWt 
addpath(genpath('../numerical methods'))


%% System definition

% Default parameter vector (B: first and second moment stable)
par = [-3.0 2.1 0 -1.5 0.0 1 1]; % [a b c alpha beta gamma tau]
tau = par(7);

% stochastic system definition
sysDDE.f = @(x,xd,t,p) p(1)*x(1,:) + p(2)*xd(1,:,1) + p(3); % drift term
sysDDE.g = @(x,xd,t,p) p(4)*x(1,:) + p(5)*xd(1,:,1) + p(6); % diffustion term
sysDDE.dg = @(x,xd,t,p) p(4); % diffustion term Jacobian

% Linear SDDE definition (with dummy delays)
sysLDDE.F = @(t,p) p(1);                % present drift term
sysLDDE.Fd = {@(t,p) p(2), @(t,p) 0};   % delayed drift term
sysLDDE.f = @(t,p) p(3);                % constant drift term
sysLDDE.G = @(t,p) p(4);                % present diffusion term
sysLDDE.Gd = {@(t,p) p(5), @(t,p) 0};   % delayed diffusion term
sysLDDE.g = @(t,p) p(6);                % constant diffusion term

% History function
yhist = @(t) 1 + 0.0*t;
M = 50; % discretization dimension
theta = cheb_mesh(M,[-3*tau,0]);

% stochastic ODE definition (pseudo-spectral approximation)
sysODE.f = @(x,t,p) sdde_ps_approx(x,t,p,sysDDE.f,1,M,tau); % extended drift term
sysODE.g = @(x,t,p) [zeros(M-1,1); sysDDE.g(x(end),x(1),t,p)]; % extended diffustion term

% Simulation options
sim.T = 50;         % simulation end time
sim.dt = 1e-3;      % timestep (increment size)
sim.n = 5000;       % number of generated realizations
sim.logs = true;    % log solver progress
N = 5000; % gather multiple samples to average the found correlations

%% Find the stationary distribution

% solution mean: (a+b)*x - c = 0
x_mean = -sysLDDE.f(0,par) / (sysLDDE.F(0,par) + sysLDDE.Fd{1}(0,par));
xt_mean = @(t) x_mean*ones(M,1);

% second moment: dC/dt = 0, JC*bC = 0
bC = pde_cov_const(sim.T,xt_mean,sysLDDE,par,[tau 3*tau],M);
JC = pde_cov_jac(0,sysLDDE,par,[tau 3*tau],M);
C_mean = reshape(-JC\reshape(bC,[],1),M,M) - xt_mean(0) * xt_mean(0).';
C_stat = C_mean(M,:);


%% Evaluate second moment based on MC simulation
C_sim = zeros(M,N);
% figure(); plot(theta,C_stat,'ko-'); hold on;
for i = 1:N
    % Run simulation using the DDE form adn evaluate the covariance
    [tu,U,Ud,~] = sdde_EM_sim(sysDDE,par,yhist,tau,sim);
    is = find(tu>=sim.T-3*tau); tui = tu(is);
    Ui = squeeze(U(1,is,:)); Vtu = cov(Ui.');
    C_sim(:,i) = interp1(tui,Vtu(end,:),theta+sim.T,'linear');
    % plot(theta,C_sim(:,i)); drawnow
    fprintf('Run %i/%i completed\n',i,N)
end
% xlabel('$\vartheta$'); ylabel('$\bar{\phi}(\vartheta)$'); hold off

%% Save simulation data
if ~exist('data','dir')
    mkdir('data')
end
save('data/scalar_corr_B.mat','theta','C_stat','C_sim')