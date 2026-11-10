clear; close all; clc;
%% Simulate trajectories of a stochastic turning process
% dx(t) = (a*x(t) + b*x(t-tau))*dt + (alpha*x(t)+beta*x(t-tau)+gamma)*dW(t)
% a = [0 1; -1-w -2*zeta]; b = [0 0; w 0]; alpha = [0 0; -sigma 0]; 
% beta = [0 0; sigma 0]; gamma = [0; sigma]; tau = 2*pi/omega
addpath(genpath('numerical methods'))


%% System definition

% System matrices
A = @(p) [0 1; -1-p(2) -2*p(1)];
B = @(p) [0 0; p(2) 0];
alpha = @(p) [0 0; -p(4)*p(2) 0];
beta = @(p) [0 0; p(4)*p(2) 0];
gamma = @(p) [0; p(4)*p(2)];

% Default parameter vector (A: first and second-moment stable)
par = [0.05 0.1 1.0 1.0]; % [zeta w omega sigma]
tau = 2*pi/par(3);

% stochastic system definition
sysDDE.f = @(x,xd,t,p) A(p)*x + B(p)*squeeze(xd(:,:,1));  % drift term
sysDDE.g = @(x,xd,t,p) alpha(p)*x + beta(p)*squeeze(xd(:,:,1)) + gamma(p); % diffustion term

% Linear SDDE definition
sysLDDE.F = @(t,p) A(p);            % present drift term
sysLDDE.Fd{1} = @(t,p) B(p);        % delayed drift term
sysLDDE.f = @(t,p) [0; 0];          % constant drift term
sysLDDE.G = @(t,p) alpha(p);        % present diffusion term
sysLDDE.Gd = @(t,p) beta(p);        % delayed diffusion term
sysLDDE.g = @(t,p) gamma(p);        % constant diffusion term

% History function
yhist = @(t) [1; 0.0] + 0*t;

% stochastic ODE definition (pseudo-spectral approximation)
sysODE.f = @(x,t,p) sdde_ps_approx(x,t,p,sysDDE.f,d,M1,tau); % extended drift term
sysODE.g = @(x,t,p) [zeros(M1-1,1); sysDDE.g(x(end),x(1),t,p)]; % extended diffustion term

% Simulation options
sim.T = 200;        % simulation end time
sim.dt = 1e-2;      % timestep (increment size)
sim.n = 1000;       % number of generated realizations
sim.logs = true;    % log solver progress


%% Evaluate second moment based on MC simulation

% Run simulation and evaluate simulated moments
[tu,U,Ud,~] = sdde_EM_sim(sysDDE,par,yhist,tau,sim);
Exu = 0*tu; Sxu = 0*tu;
for i = 1:length(tu)
    pxU_temp = fitdist(squeeze(U(1,i,:)),'normal');
    Exu(i) = pxU_temp.mu; Sxu(i) = pxU_temp.sigma;
end
Uxp = squeeze(U(1,:,1:20));


%% 1D PDE solution of the correlation function
% 
% % deterministic solution
% M = 10; % Chebyshev mesh dimension
% tM = cheb_mesh(M,[-tau 0]);
% X0 = reshape(yhist(tM),[],1);
% sysODE = sys_ps_approx_ode(sysLDDE,tau,2,M);
% solO = ode45(@(t,y) sysODE.f(y,t,par),[0 sim.T],X0);
% xref = @(t) deval(solO,t);
% 
% % PDE time integration
% Phi0 = @(t) reshape(yhist(t)*reshape(yhist(t+tM),1,[]),[],1);
% pfun = @(t,y,Z) pde_corr(t,y,Z,xref,sysLDDE,par,tau);
% solP = dde23(pfun,-tM(1:end-1),Phi0,[0 sim.T]);
% tP = solP.x; XP = xref(tP); 
% ExP = XP(end-1,:); SxP = sqrt(abs(solP.y(end-3,:)-ExP.^2));


%% Compare simulation results

figure('Position',[400 400 400 400]); ax1 = axes(); box on; hold on
plot(tu,Uxp,'Color',[0 0 1 0.3],'HandleVisibility','off');
% plot(tP,ExP,'r',tP,ExP+SxP,'k',tP,ExP-SxP,'k','LineWidth',1.5)
plot(tu,Exu,'r',tu,Exu+Sxu,'k',tu,Exu-Sxu,'k','LineWidth',1.5)
xlabel('$t$'); ylabel('$x(t)$'); hold off
% legend('$\mu(t)$','$\mu(t)\pm\sigma(t)$');
title(sprintf('$\\tilde{w}$ = %0.2f, $\\tilde{\\omega}$ = %0.2f',...
    par(2),par(3)))
ax1.Units = 'normalized'; pos = ax1.Position; 
annotation(gcf, 'textbox', [0, pos(2)+pos(4), 0.04, 0.04], ...
    'String', '(A)','LineStyle', 'none','Interpreter','latex');
fontsize(16,"points"); grid on;
% exportgraphics(gcf,'plots/turning_sim_a.png','Resolution',400)