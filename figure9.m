clear; close all; clc;
%% Simulate the covariance of a linear, scalar SDDE
% dx(t) = (a*x(t)+b*x(t-tau) + c)*dt + 
%   (alpha*x(t)+beta*x(t-tau) + gamma)*dWt 
addpath(genpath('numerical methods'))

% RUN "scalar_sim_batch_A.m" BEFORE THIS SCIPT!


%% System definition

% Default parameter vector (A: first and second-moment stable)
par = [-3.0 1.0 0 -1.5 0.0 1 1]; % [a b c alpha beta gamma tau]
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
sim.T = 20;         % simulation end time
sim.dt = 1e-3;      % timestep (increment size)
sim.n = 5000;       % number of generated realizations
sim.logs = true;    % log solver progress
n_step = 100;       % postprocessing timesteps


%% Evaluate second moment based on MC simulation

% Run simulation using the DDE form adn evaluate the covariance
[tu,U,Ud,~] = sdde_EM_sim(sysDDE,par,yhist,tau,sim);
% [tu,U,Ud,~] = sdde_MI_sim(sysDDE,par,yhist,tau,sim);
% UT = [squeeze(Ud(1,end,:,:)) squeeze(U(1,end,:))];
% CU = cov(UT);
% figure(); surf(theta,theta,CU)

% Evaluate quartiles and mean
UP = prctile(squeeze(U), 25, 2);
U0 = prctile(squeeze(U), 50, 2);
UM = prctile(squeeze(U), 75, 2);
UA = mean(squeeze(U),2);

% Evaluate mean and standard deviation
% Eu = mean(squeeze(U),2);
% Su = std(squeeze(U),[],2);
% UP = Eu+Su; U0 = Eu; UM = Eu-Su;

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

% Final extended covariance matrix
CUTe = cov(squeeze(U(1,end-3/sim.dt:end,:)).');
thCe = linspace(-3,0,size(CUTe,1));


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


%% Evaluate first and second moment stability

% first moment
sysLODE = sys_ps_approx_lode(sysLDDE,tau,1,M);
A = sysLODE.F(0,par);

% second moment
JC = pde_cov_jac(0,sysLDDE,par,tau,M);
% dC = JC*solC.y(:,end)-pde_cov(sim.T,solC.y(:,end),xref,sysLDDE,par,tau);
% disp(max(abs(dC),[],"all"))
% figure(); plot(dC);


%% Find the stationary distribution

% solution mean: (a+b)*x - c = 0
x_mean = -sysLDDE.f(0,par) / (sysLDDE.F(0,par) + sysLDDE.Fd(0,par));
xt_mean = @(t) x_mean*ones(M,1);

% second moment: dC/dt = 0, JC*bC = 0
bC = pde_cov_const(sim.T,xt_mean,sysLDDE,par,tau,M);
C_mean = reshape(-JC\reshape(bC,[],1),M,M) - xt_mean(0) * xt_mean(0).';
s_mean = sqrt(C_mean(end,end));


%% Analytic solution of the stationary distribution
% For now, restricted to c = 0, x_mean = 0!

% BVP solution on the [-\tau 0] interval
phi_a = corr_scalar_asol(sysLDDE,tau,par);
[th_m, phi_m] = corr_stat_asol(sysLDDE,tau,par);
% figure(); plot(th_m,phi_a(th_m),th_m,phi_m,':')

% DDE IVP solution on the theta<-\tau domain
y0 = @(t) phi_a(-t-tau(1)); tau_max= 3;
sol = dde23(@(t,y,Z)corr_stat_dde(t,y,Z,sysLDDE,par,x_mean),tau,...
    y0,[0 tau_max-tau(1)]);
phi_n = @(th) deval(sol,-th-tau(1));

% Combine the extended analytic solution
Mt = 1000; Ca_ext = zeros(Mt); 
theta2 = linspace(-tau_max,0,Mt); 
for i = 1:Mt
    for j = 1:Mt
        th_ij = -abs(theta2(i) - theta2(j));
        if th_ij < -tau(1)
            Ca_ext(i,j) = phi_n(th_ij);
        else
            Ca_ext(i,j) = phi_a(th_ij);
        end
    end
end
Ca_ext = Ca_ext - x_mean^2;


%% Combined plot of results
figure('Position',[400 400 950 500]);

subplot(2,3,1); box on; hold on; grid on;
plot(tu,UA,'r','LineWidth',1.5);
plot(tu,U0,'b','LineWidth',1.5);
plot(tu,UP,'m','LineWidth',1.5);
plot(tu,UM,'m','LineWidth',1.5,'HandleVisibility','off');
yline(0,'k:','LineWidth',1.5,'HandleVisibility','off')
% fill([tu, fliplr(tu)], [UM.', fliplr(UP')], [0 0 1], ...
%     'EdgeColor', 'none', 'FaceAlpha', 0.3);
% title(sprintf('$a$ = %0.1f, $b$ = %0.1f',par(1),par(2)))
xlabel('$t$'); ylabel('$x$')
annotation(gcf, 'textbox', [0.07, 0.92, 0.05, 0.05], ...
    'String', '(a)','LineStyle', 'none','Interpreter','latex');
% legend('mean($x$)','$q_{0.5}(x)$','$q_{0.25,0.75}(x)$')

% subplot(2,3,1); box on; hold on; grid on;
% plot(tu,UA,'b','LineWidth',1.5);
% plot(tu,U0,'k','LineWidth',1.5);
% p = patch([tu, fliplr(tu)],[UP.', fliplr(UM.')], [0 0 1]);  
% p.FaceAlpha = 0.2; p.EdgeColor = 'none';
% yline(0,'k:','LineWidth',1.5,'HandleVisibility','off')
% % fill([tu, fliplr(tu)], [UM.', fliplr(UP')], [0 0 1], ...
% %     'EdgeColor', 'none', 'FaceAlpha', 0.3);
% % title(sprintf('$a$ = %0.1f, $b$ = %0.1f',par(1),par(2)))
% xlabel('$t$'); ylabel('$x$')
% annotation(gcf, 'textbox', [0.07, 0.92, 0.05, 0.05], ...
%     'String', '(a)','LineStyle', 'none','Interpreter','latex');
% % legend('mean($x$)','$q_{0.5}(x)$','$q_{0.25,0.75}(x)$')

subplot(2,3,2); plot_lambda(eig(A)); grid on;
s = findobj(gca, 'Type', 'Scatter'); 
xlabel('Re($\lambda$)'); ylabel('Im($\lambda$)')
s(1).LineWidth = 1.5; s(2).LineWidth = 1.5;
s(1).Marker = '.'; s(2).Marker = '.';
s(1).SizeData = 75; s(2).SizeData = 75;
annotation(gcf, 'textbox', [0.35, 0.92, 0.05, 0.05], ...
    'String', '(b)','LineStyle', 'none','Interpreter','latex');
% title('');

subplot(2,3,3); plot_lambda(eig(JC));  grid on;
xlabel('Re($\lambda$)'); ylabel('Im($\lambda$)')
s = findobj(gca, 'Type', 'Scatter'); 
s(1).LineWidth = 1.5; s(2).LineWidth = 1.5;
s(1).Marker = '.'; s(2).Marker = '.';
s(1).SizeData = 75; s(2).SizeData = 75;
annotation(gcf, 'textbox', [0.63, 0.92, 0.05, 0.05], ...
    'String', '(c)','LineStyle', 'none','Interpreter','latex');
% title('');

subplot(2,3,4); box on; hold on; grid on;
plot(tus,Vtu,'b',tc,Ct(end,:),'r','LineWidth',1.5);
yline(s_mean^2,'k:','LineWidth',1.5)
xlabel('$t$'); ylabel('$V[x]$')
annotation(gcf, 'textbox', [0.07, 0.43, 0.05, 0.05], ...
    'String', '(d)','LineStyle', 'none','Interpreter','latex');

subplot(2,3,5); box on; hold on;  grid on;
% surf(theta,theta,CUT,'EdgeColor','b','EdgeAlpha',0.1, ...
%     'FaceColor','b','FaceAlpha',0.2);
% surf(thCe,thCe,CUTe,'EdgeColor','none', ...
%     'FaceColor','b','FaceAlpha',0.2);
% surf(theta,theta,C_mean,'EdgeColor','none', ...
%     'FaceColor','r','FaceAlpha',0.9);
surf(theta2,theta2,Ca_ext,'EdgeColor','none', ...
    'FaceColor','r','FaceAlpha',0.9);
xlabel('$\theta$'); ylabel('$\vartheta$'); 
zlabel('$\bar{C}(\theta,\vartheta)$'); hold off
view(-30,30)
lightangle(45,45); light('Position',[-1 1 1]);
lighting gouraud; material shiny
annotation(gcf, 'textbox', [0.35, 0.43, 0.05, 0.05], ...
    'String', '(e)','LineStyle', 'none','Interpreter','latex');

% subplot(2,3,6); box on; hold on;  grid on;
% plot(theta2,Ca_ext(end,:),'r-', thCe,CUTe(end,:),'b-',...
%     'LineWidth',1.5,'MarkerSize',10);
% xlabel('$\vartheta$'); ylabel('$\bar{\phi}(\vartheta)$'); 
% xline(-tau,':k','HandleVisibility','off'); xlim('tight'); hold off
% annotation(gcf, 'textbox', [0.63, 0.43, 0.05, 0.05], ...
%     'String', '(f)','LineStyle', 'none','Interpreter','latex');

load('data/scalar_corr_A.mat','theta','C_stat','C_sim')
mC = mean(C_sim,2).'; sC = sqrt(var(C_sim,[],2)).';
subplot(2,3,6); box on; hold on;  grid on;
p = patch([theta, fliplr(theta)],[mC-sC, fliplr(mC+sC)], [0 0 1]);  
p.FaceAlpha = 0.2; p.EdgeColor = 'none';
plot(theta,C_stat,'r-', theta,mean(C_sim,2),'b-',...
    'LineWidth',1.5,'MarkerSize',10);
xlabel('$\vartheta$'); ylabel('$\bar{\phi}(\vartheta)$'); ylim([0 0.3]);
xline(-tau,':k','HandleVisibility','off'); xlim('tight'); hold off
annotation(gcf, 'textbox', [0.63, 0.43, 0.05, 0.05], ...
    'String', '(f)','LineStyle', 'none','Interpreter','latex');
% annotation(gcf, 'textbox', [0.05, 0.9, 0.05, 0.05], ...
%     'String', '(A)','LineStyle', 'none','Interpreter','latex');
fontsize(16,"points")
exportgraphics(gcf,'plots/scalar_sim_a.png','Resolution',400)