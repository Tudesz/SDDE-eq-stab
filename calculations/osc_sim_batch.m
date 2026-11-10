clear; close all; clc;
%% Simulate the covariance of a stochastic 2D linear equation
% dx(t) = (A*x(t) + B*x(t-tau))*dt + (alpha*x(t)+beta*x(t-tau)+gamma)*dW(t)
% fodor2023: A = [0 1; k 0]; B = [0 0; -p -d];
%   alpha = [0 0; 0 0]; beta = -sigma*B; gamma = [0; sigma]
addpath(genpath('../numerical methods'))


%% System definition

% System matrices (fodor2023)
A = @(p) [0 1; p(4) 0];
B = @(p) [0 0; -p(1) -p(2)];
alpha = @(p) zeros(2);
beta = @(p) [0 0; p(3)*p(1) p(3)*p(2)];
gamma = @(p) [0; p(3)];

% Default parameter vector
par = [6.5 3.5 0.1 5 0.3]; % [p d sigma k tau]
tau_max = 1; tau = [par(5) tau_max]; d = 2;

% stochastic system definition
sysDDE.f = @(x,xd,t,p) A(p)*x + B(p)*squeeze(xd(:,:,1));  % drift term
sysDDE.g = @(x,xd,t,p) alpha(p)*x + beta(p)*squeeze(xd(:,:,1)) + gamma(p); % diffustion term

% Linear SDDE definition
sysLDDE.F = @(t,p) A(p);            % present drift term
sysLDDE.Fd{1} = @(t,p) B(p);        % delayed drift term
sysLDDE.Fd{2} = @(t,p) zeros(2);    % delayed drift term
sysLDDE.f = @(t,p) [0; 0];          % constant drift term
sysLDDE.G = @(t,p) alpha(p);        % present diffusion term
sysLDDE.Gd{1} = @(t,p) beta(p);     % delayed diffusion term
sysLDDE.Gd{2} = @(t,p) zeros(2);    % present diffusion term
sysLDDE.g = @(t,p) gamma(p);        % constant diffusion term

% History function
yhist = @(t) [1.0; 0.0] + 0*t;
M = 50; % simulation output discretization
theta = cheb_mesh(M,[-tau_max,0]);

% Simulation options
sim.T = 10;         % simulation end time
sim.dt = 1e-3;      % timestep (increment size)
sim.n = 1000;       % number of generated realizations
sim.logs = true;    % log solver progress
N = 1000; % gather multiple samples to average the found correlations


%% Find the stationary distribution

% solution mean: (a+b)*x - c = 0
x_mean = -(sysLDDE.F(0,par) + sysLDDE.Fd{1}(0,par))\sysLDDE.f(0,par);
xt_mean = @(t) kron(ones(M,1),x_mean);

% second moment: dC/dt = 0, JC*C = bC
JC = pde_cov_jac(0,sysLDDE,par,tau,M);
bC = pde_cov_const(0,xt_mean,sysLDDE,par,tau,M);
C_mean = reshape(-JC\reshape(bC,[],1),d*M,d*M) - xt_mean(0) * xt_mean(0).';
C_stat = [C_mean(1:2:end,end-1) C_mean(2:2:end,end-1)...
    C_mean(1:2:end,end) C_mean(2:2:end,end)];

%% Evaluate second moment based on MC simulation
C_sim = zeros(4,M,N);
% figure(); plot(theta,C_stat,'ko-'); hold on;
for i = 1:N
    % Run simulation using the DDE form adn evaluate the covariance
    [tu,U,Ud,~] = sdde_EM_sim(sysDDE,par,yhist,tau,sim);
    is = find(tu>=sim.T-tau_max); tui = tu(is);
    Ui = reshape(U(:,is,:),[],sim.n); 
    Vtu = cov(Ui.');
    C_sim(1,:,i) = interp1(tui,Vtu(end-1,1:2:end),theta+sim.T,'linear');
    C_sim(2,:,i) = interp1(tui,Vtu(end-1,2:2:end),theta+sim.T,'linear');
    C_sim(3,:,i) = interp1(tui,Vtu(end,1:2:end),theta+sim.T,'linear');
    C_sim(4,:,i) = interp1(tui,Vtu(end,2:2:end),theta+sim.T,'linear');
    % plot(theta,squeeze(C_sim(:,:,i)),'b'); drawnow
    fprintf('Run %i/%i completed\n',i,N)
end
% xlabel('$\vartheta$'); ylabel('$\bar{\phi}(\vartheta)$'); hold off

%% Save simulation data
save('../data/osc_corr.mat','theta','C_stat','C_sim')

