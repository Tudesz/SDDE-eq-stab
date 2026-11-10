clear; close all; clc;
%% Simulate the covariance of a stochastic 2D linear equation
% dx(t) = (A*x(t) + B*x(t-tau))*dt + (alpha*x(t)+beta*x(t-tau)+gamma)*dW(t)
% fodor2023: A = [0 1; k 0]; B = [0 0; -p -d];
%   alpha = [0 0; 0 0]; beta = -sigma*B; gamma = [0; sigma]
addpath(genpath('numerical methods'))

% RUN "osc_sim_batch.m" BEFORE THIS SCIPT!

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
M1 = 40; % simulation output discretization
theta1 = cheb_mesh(M1,[-tau_max,0]);
M2 = 20; % stationary distribution discretization
theta2 = cheb_mesh(M2,[-tau_max,0]);

% stochastic ODE definition (pseudo-spectral approximation)
sysODE.f = @(x,t,p) sdde_ps_approx(x,t,p,sysDDE.f,d,M1,tau); % extended drift term
sysODE.g = @(x,t,p) [zeros(M1-1,1); sysDDE.g(x(end),x(1),t,p)]; % extended diffustion term

% Simulation options
sim.T = 10;         % simulation end time
sim.dt = 1e-3;      % timestep (increment size)
sim.n = 1000;       % number of generated realizations
sim.logs = true;    % log solver progress


%% Evaluate second moment based on MC simulation

% Run simulation using the DDE form adn evaluate the covariance
[tu,U,Ud,~] = sdde_EM_sim(sysDDE,par,yhist,[tau(1) -theta1(1:end-1)],sim);
UT = squeeze(cat(1,reshape(permute(Ud(:,end,:,2:end),[1 4 2 3]),...
    [],1,size(U,3)),U(:,end,:))).';
CU = cov(UT);
% figure(); plot(tu,squeeze(U(1,:,1:20)),'b')
% figure(); surf(theta1,theta1,CU(2:d:end,2:d:end))
Exu = 0*tu; Sxu = 0*tu;
Evu = 0*tu; Svu = 0*tu;
for i = 1:length(tu)
    pxU_temp = fitdist(squeeze(U(1,i,:)),'normal');
    Exu(i) = pxU_temp.mu; Sxu(i) = pxU_temp.sigma;
    pvU_temp = fitdist(squeeze(U(2,i,:)),'normal');
    Evu(i) = pvU_temp.mu; Svu(i) = pvU_temp.sigma;
end


%% Discretized stationary distribution

% solution mean: (a+b)*x - c = 0
x_mean = -(sysLDDE.F(0,par) + sysLDDE.Fd{1}(0,par))\sysLDDE.f(0,par);
xt_mean = @(t) kron(ones(M2,1),x_mean);

% second moment: dC/dt = 0, JC*C = bC
% JC = pde_cov_jac(0,sysLDDE,par,tau,M2);
% bC = pde_cov_const(0,xt_mean,sysLDDE,par,tau,M2);
% C_mean = reshape(-JC\reshape(bC,[],1),d*M2,d*M2) - xt_mean(0) * xt_mean(0).';

% correaltaiotn function: dPhi/dt = 0, JPh*Ph = bPh
% JPh = pde_corr_jac(0,sysLDDE,par,tau,M2);
% bPh = pde_corr_const(0,xt_mean,sysLDDE,par,tau,M2);
% Ph_mean = reshape(-JPh\bPh,d^2,M2);
% figure(); plot(theta2,C_mean(end,2:d:end),theta2,Ph_mean(end,:),':')


%% Analytic solution of the stationary distribution
% For now, restricted to c = 0, x_mean = 0!

% BVP solution on the [-tau 0] interval
[th_a,phi_a] = corr_stat_asol(sysLDDE,tau,par);
phi_i = @(th) interp1(th_a,phi_a.',th,'spline').';
% figure(); plot(th_a,phi_a)

% DDE IVP solution on the theta<-\tau domain
y0 = @(t) phi_i(-t-tau(1));
if tau_max > tau(1)
    sol = dde23(@(t,y,Z)corr_stat_dde(t,y,Z,sysLDDE,par,x_mean),tau,...
        y0,[0 tau_max-tau(1)]);
end
phi_n = @(th) deval(sol,-th-tau(1));

% Combine the extended analytic solution
Mt = 100; Ca_ext = zeros(d*Mt);
theta3 = linspace(-tau_max,0,Mt); 
for i = 1:Mt
    ii = (i-1)*d+1:i*d;
    for j = 1:Mt
        jj = (j-1)*d+1:j*d;
        th_ij = -abs(theta3(i) - theta3(j));
        if th_ij < -tau(1)
            Ca_ext(ii,jj) = reshape(phi_n(th_ij),d,d);
        else
            Ca_ext(ii,jj) = reshape(phi_i(th_ij),d,d);
        end
    end
end
sv_mean = sqrt(Ca_ext(end,end));
sx_mean = sqrt(Ca_ext(end-1,end-1));


%% Combined plot of results

% figure('Position',[400 400 800 400]);
figure('Position',[400 400 600 500]);

% ax1 = subplot(2,2,1); box on; hold on
ax1 = subplot(2,1,1); box on; hold on;
plot(tu,squeeze(U(1,:,1:20)),'Color',[0 0 1 0.2],'HandleVisibility','off')
plot(tu,Exu,'r',tu,Exu+Sxu,'k','LineWidth',1.5);
plot(tu,Exu-Sxu,'k','LineWidth',1.5,'HandleVisibility','off');
yline(x_mean(1),'--k','LineWidth',1)
yline(x_mean(1)+[-sx_mean sx_mean],':k','LineWidth',1)
xlabel('$t$'); ylabel('$x(t)$'); hold off; grid on;
ylim([-0.35 1.1])
% legend('$\mu(t)$','$\mu(t)\pm\sigma(t)$');
ax1.Units = 'normalized'; pos = ax1.Position; 
% annotation(gcf, 'textbox', [pos(1)-0.07, pos(2)+pos(4), 0.04, 0.04], ...
%     'String', '(a)','LineStyle', 'none','Interpreter','latex');
annotation(gcf, 'textbox', [0, pos(2)+pos(4), 0.05, 0.05], ...
    'String', '(a)','LineStyle', 'none','Interpreter','latex');

% ax2 = subplot(2,2,3); box on; hold on
ax2 = subplot(2,1,2); box on; hold on
plot(tu,squeeze(U(2,:,1:20)),'Color',[0 0 1 0.1],'HandleVisibility','off')
plot(tu,Evu,'r',tu,Evu+Svu,'k','LineWidth',1.5);
plot(tu,Evu-Svu,'k','LineWidth',1.5,'HandleVisibility','off');
yline(x_mean(2),'--k','LineWidth',1)
yline(x_mean(2)+[-sv_mean sv_mean],':k','LineWidth',1)
xlabel('$t$'); ylabel('$v(t)$'); hold off; grid on;
ylim([-1.5 1.0])
ax2.Units = 'normalized'; pos = ax2.Position; 
% annotation(gcf, 'textbox', [pos(1)-0.07, pos(2)+pos(4), 0.04, 0.04], ...
%     'String', '(b)','LineStyle', 'none','Interpreter','latex');
annotation(gcf, 'textbox', [0, pos(2)+pos(4), 0.05, 0.05], ...
    'String', '(b)','LineStyle', 'none','Interpreter','latex');
fontsize(20,"points")
% exportgraphics(gcf,'plots/osc_stationary1.png','Resolution',400)

% % ax3 = axes('Position',[0.62 0.17 0.3347 0.75]);  box on; hold on
% figure('Position',[400 400 600 500]); ax3 = axes();  box on; hold on
% x_CC = C_mean(end-1,1:2:end);
% x_MC = CU(end-1,1:2:end);
% x_A = Ca_ext(end-1,1:2:end);
% plot(theta2,x_CC,'k-','HandleVisibility','on','LineWidth',1.5,'MarkerSize',10);
% plot(theta1,x_MC,'b-','HandleVisibility','on','LineWidth',1.5,'MarkerSize',10);
% plot(theta3,x_A,'r-','HandleVisibility','on','LineWidth',1.5,'MarkerSize',10);
% v_CC = C_mean(end,2:2:end);
% v_MC = CU(end,2:2:end);
% v_A = Ca_ext(end,2:2:end);
% plot(theta2,v_CC,'k--','HandleVisibility','off','LineWidth',1.5,'MarkerSize',10);
% plot(theta1,v_MC,'b--','HandleVisibility','off','LineWidth',1.5,'MarkerSize',10);
% plot(theta3,v_A,'r--','HandleVisibility','off','LineWidth',1.5,'MarkerSize',10);
% xv_CC = C_mean(end,1:2:end);
% xv_MC = CU(end,1:2:end);
% xv_A = Ca_ext(end,1:2:end);
% plot(theta2,xv_CC,'k:','HandleVisibility','off','LineWidth',1.5,'MarkerSize',10);
% plot(theta1,xv_MC,'b:','HandleVisibility','off','LineWidth',1.5,'MarkerSize',10);
% plot(theta3,xv_A,'r:','HandleVisibility','off','LineWidth',1.5,'MarkerSize',10);
% vx_CC = C_mean(end-1,2:2:end);
% vx_MC = CU(end-1,2:2:end);
% vx_A = Ca_ext(end-1,2:2:end);
% plot(theta2,vx_CC,'k-.','HandleVisibility','off','LineWidth',1.5,'MarkerSize',10);
% plot(theta1,vx_MC,'b-.','HandleVisibility','off','LineWidth',1.5,'MarkerSize',10);
% plot(theta3,vx_A,'r-.','HandleVisibility','off','LineWidth',1.5,'MarkerSize',10);
% legend('Discretization','MC simulation','Analytic','Location','northwest')
% xlabel('$\vartheta$'); ylabel('$\bar{\phi}(\vartheta)$'); 
% xline(-tau,':k','HandleVisibility','off'); xlim('tight'); hold off
% ax3.Units = 'normalized'; pos = ax3.Position; 
% % annotation(gcf, 'textbox', [pos(1)-0.12, pos(2)+pos(4), 0.04, 0.04], ...
% %     'String', '(c)','LineStyle', 'none','Interpreter','latex');
% annotation(gcf, 'textbox', [0, pos(2)+pos(4), 0.05, 0.05], ...
%     'String', '(c)','LineStyle', 'none','Interpreter','latex');
% fontsize(16,"points")
% exportgraphics(gcf,'plots/osc_stationary2.png','Resolution',400)

%%
load('data\osc_corr.mat','theta','C_stat','C_sim')
mC = mean(C_sim,3).'; sC = sqrt(var(C_sim,[],3)).';
% ax3 = axes('Position',[0.62 0.17 0.3347 0.75]);  box on; hold on
figure('Position',[400 400 600 500]); ax3 = axes();  box on; hold on
x_CC = C_stat(:,1); x_MC = mC(:,1); x_A = Ca_ext(end-1,1:2:end);
plot(theta,x_CC,'k-','HandleVisibility','on','LineWidth',1.5,'MarkerSize',10);
plot(theta,x_MC,'b-','HandleVisibility','on','LineWidth',1.5,'MarkerSize',10);
plot(theta3,x_A,'r-','HandleVisibility','on','LineWidth',1.5,'MarkerSize',10);
p = patch([theta, fliplr(theta)],[mC(:,1).'-sC(:,1).', fliplr(mC(:,1).'+sC(:,1).')], [0 0 1]);  
p.FaceAlpha = 0.2; p.EdgeColor = 'none'; p.HandleVisibility = 'off';
v_CC = C_stat(:,4); v_MC = mC(:,4); v_A = Ca_ext(end,2:2:end);
plot(theta,v_CC,'k--','HandleVisibility','off','LineWidth',1.5,'MarkerSize',10);
plot(theta,v_MC,'b--','HandleVisibility','off','LineWidth',1.5,'MarkerSize',10);
plot(theta3,v_A,'r--','HandleVisibility','off','LineWidth',1.5,'MarkerSize',10);
p = patch([theta, fliplr(theta)],[mC(:,4).'-sC(:,4).', fliplr(mC(:,4).'+sC(:,4).')], [0 0 1]);  
p.FaceAlpha = 0.2; p.EdgeColor = 'none'; p.HandleVisibility = 'off';
xv_CC = C_stat(:,3); xv_MC = mC(:,3); xv_A = Ca_ext(end,1:2:end);
plot(theta,xv_CC,'k:','HandleVisibility','off','LineWidth',1.5,'MarkerSize',10);
plot(theta,xv_MC,'b:','HandleVisibility','off','LineWidth',1.5,'MarkerSize',10);
plot(theta3,xv_A,'r:','HandleVisibility','off','LineWidth',1.5,'MarkerSize',10);
p = patch([theta, fliplr(theta)],[mC(:,3).'-sC(:,3).', fliplr(mC(:,3).'+sC(:,3).')], [0 0 1]);  
p.FaceAlpha = 0.2; p.EdgeColor = 'none'; p.HandleVisibility = 'off';
vx_CC = C_stat(:,2); vx_MC = mC(:,2); vx_A = Ca_ext(end-1,2:2:end);
plot(theta,vx_CC,'k-.','HandleVisibility','off','LineWidth',1.5,'MarkerSize',10);
plot(theta,vx_MC,'b-.','HandleVisibility','off','LineWidth',1.5,'MarkerSize',10);
plot(theta3,vx_A,'r-.','HandleVisibility','off','LineWidth',1.5,'MarkerSize',10);
legend('Discretization','MC simulation','Analytic','Location','northwest')
p = patch([theta, fliplr(theta)],[mC(:,2).'-sC(:,2).', fliplr(mC(:,2).'+sC(:,2).')], [0 0 1]);  
p.FaceAlpha = 0.2; p.EdgeColor = 'none'; p.HandleVisibility = 'off';
xlabel('$\vartheta$'); ylabel('$\bar{\phi}(\vartheta)$'); 
xline(-tau,':k','HandleVisibility','off'); xlim('tight'); hold off; grid on;
ax3.Units = 'normalized'; pos = ax3.Position; 
% annotation(gcf, 'textbox', [pos(1)-0.12, pos(2)+pos(4), 0.04, 0.04], ...
%     'String', '(c)','LineStyle', 'none','Interpreter','latex');
annotation(gcf, 'textbox', [0, pos(2)+pos(4), 0.05, 0.05], ...
    'String', '(c)','LineStyle', 'none','Interpreter','latex');
fontsize(20,"points")
% exportgraphics(gcf,'plots/osc_stationary2.png','Resolution',400)

% exportgraphics(gcf,'osc_stationary.png','Resolution',400)

