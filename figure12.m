clear; close all; clc;
%% Simulate the covariance of a linear, scalar SDDE
% dx(t) = (a*x(t)+b*x(t-tau) + c)*dt + 
%   (alpha*x(t)+beta*x(t-tau) + gamma)*dWt 
addpath(genpath('numerical methods'))


%% System definition

% Default parameter vector (D: first and second-moment unstable)
par = [1 -1.5 0 -1.5 0.0 1 1]; % [a b c alpha beta gamma tau]
tau = par(7);

% stochastic system definition
sysDDE.f = @(x,xd,t,p) p(1)*x(1,:) + p(2)*xd(1,:,1) + p(3); % drift term
sysDDE.g = @(x,xd,t,p) p(4)*x(1,:) + p(5)*xd(1,:,1) + p(6); % diffustion term
sysDDE.dg = @(x,xd,t,p) p(4); % diffustion term Jacobian

% Linear SDDE definition
sysLDDE.F = @(t,p) p(1);    % present drift term
sysLDDE.Fd = @(t,p) p(2);   % delayed drift term
sysLDDE.f = @(t,p) p(3);    % constant drift term
sysLDDE.G = @(t,p) p(4);    % present diffusion term
sysLDDE.Gd = @(t,p) p(5);   % delayed diffusion term
sysLDDE.g = @(t,p) p(6);    % constant diffusion term

% History function
yhist = @(t) 1 + 0.0*t;
M = 20; % discretization dimension
theta = cheb_mesh(M,[-tau,0]);

% stochastic ODE definition (pseudo-spectral approximation)
sysODE.f = @(x,t,p) sdde_ps_approx(x,t,p,sysDDE.f,1,M,tau); % extended drift term
sysODE.g = @(x,t,p) [zeros(M-1,1); sysDDE.g(x(end),x(1),t,p)]; % extended diffustion term

% Simulation options
sim.T = 2.5;        % simulation end time
sim.dt = 1e-3;      % timestep (increment size)
sim.n = 5000;       % number of generated realizations
sim.logs = true;    % log solver progress
n_step = 100;       % postprocessing timesteps


%% Evaluate second moment based on MC simulation

% Run simulation using the DDE form adn evaluate the covariance
% [tu,U,Ud,~] = sdde_EM_sim(sysDDE,par,yhist,-theta(1:end-1),sim);
[tu,U,Ud,~] = sdde_MI_sim(sysDDE,par,yhist,tau,sim);
% UT = [squeeze(Ud(1,end,:,:)) squeeze(U(1,end,:))];
% CU = cov(UT);
% figure(); surf(theta,theta,CU)

% Evaluate quartiles and mean
% UP = prctile(squeeze(U), 15.87, 2);
% U0 = prctile(squeeze(U), 50, 2);
% UM = prctile(squeeze(U), 84.13, 2);

% Evaluate mean and standard deviation
% Eu = mean(squeeze(U),2);
% Su = std(squeeze(U),[],2);
% UP = Eu+Su; U0 = Eu; UM = Eu-Su;

% Evaluate quartiles and mean
UP = prctile(squeeze(U), 25, 2);
U0 = prctile(squeeze(U), 50, 2);
UM = prctile(squeeze(U), 75, 2);
UA = mean(squeeze(U),2);

% Find Covariance operators
is = 1:n_step:length(tu);
tus = tu(is); Vtu = 0*tus;
% CUts = zeros(M^2,length(is));
for i = 1:length(tus)
    UTi = [squeeze(Ud(1,is(i),:,:)) squeeze(U(1,is(i),:))];
    % CUts(:,i) = reshape(UTi.'*UTi/sim.n,[],1);
    % CUts(:,i) = reshape(cov(UTi),[],1);
    Vtu(i) = cov(squeeze(U(1,is(i),:)));
end
% figure(); plot(tus,CUts,tus,CUts(end,:),'r')
% CUT = reshape(CUts(:,end),M,M);


%% 2D PDE solution of the second moment "function"

% deterministic solution
solO = ode45(@(t,y) sysODE.f(y,t,par),[0 sim.T],yhist(theta));
xref = @(t) deval(solO,t); xrefT = xref(sim.T);
% figure(); plot(solO.x,solO.y)

% Time integration
C0 = reshape(yhist(theta).'*yhist(theta),[],1);
cfun = @(t,y) pde_cov(t,y,xref,sysLDDE,par,tau);
solC = ode45(cfun,[0 sim.T],C0);
CD = reshape(solC.y(:,end),M,M) - xrefT * xrefT.';
tc = solC.x; Ct =zeros(M^2,length(tc));
for i = 1:length(tc)
    Ct(:,i) = solC.y(:,i) - reshape(xref(tc(i))*xref(tc(i)).',[],1);
end
% figure(); plot(tc,Ct,tc,Ct(end,:),'r')


%% 1D PDE solution of the correlation function

% Time integration
Phi0 = @(t) reshape(yhist(t)*reshape(yhist(t+theta),1,[]),[],1);
pfun = @(t,y,Z) pde_corr(t,y,Z,xref,sysLDDE,par,tau);
solP = dde23(pfun,-theta(1:end-1),Phi0,[0 sim.T]);
Phi1 = deval(solP,sim.T+theta);
CP = cov_interp(reshape(Phi1.',M,M),theta) - xrefT * xrefT.';
x_end  = @(x) x(end,:);
Phi = solP.y - x_end(xref(solP.x)) .* xref(solP.x);


%% Combined plot of results
figure('Position',[400 400 630 250]);

subplot(1,2,1); box on; hold on; grid on;
plot(tu,UA,'r','LineWidth',1.5);
plot(tu,U0,'b','LineWidth',1.5);
plot(tu,UP,'m','LineWidth',1.5);
plot(tu,UM,'m','LineWidth',1.5,'HandleVisibility','off');
yline(0,'k:','LineWidth',1.5,'HandleVisibility','off')
% fill([tu, fliplr(tu)], [UM.', fliplr(UP')], [0 0 1], ...
%     'EdgeColor', 'none', 'FaceAlpha', 0.3);
% title(sprintf('$a$ = %0.1f, $b$ = %0.1f',par(1),par(2)))
xlabel('$t$'); ylabel('$x$')
annotation(gcf, 'textbox', [0.05, 0.92, 0.05, 0.05], ...
    'String', '(a)','LineStyle', 'none','Interpreter','latex');
% axis([0 sim.T -20 20])

subplot(1,2,2); box on; hold on; grid on
plot(tus,Vtu,'b',tc,Ct(end,:),'r','LineWidth',1.5);
xlabel('$t$'); ylabel('$V[x]$')
annotation(gcf, 'textbox', [0.47, 0.92, 0.05, 0.05], ...
    'String', '(b)','LineStyle', 'none','Interpreter','latex');
axis([0 sim.T 0 1000])

% annotation(gcf, 'textbox', [0.05, 0.9, 0.05, 0.05], ...
%     'String', '(D)','LineStyle', 'none','Interpreter','latex');
fontsize(16,"points")
exportgraphics(gcf,'plots/scalar_sim_d.png','Resolution',400)