clear; close all; clc;
%% Stationary solutions to a linear, scalar SDDE
% dx(t) = (a*x(t)+b*x(t-tau) + c)*dt + 
%   (alpha*x(t)+beta*x(t-tau) + gamma)*dWt 
addpath(genpath('numerical methods'))


%% System definition

% Default parameter vector
par = [0 -2 0 -1.5 0.5 1 1]; % [a b c alpha beta gamma tau]
tau = par(7); d = 1;

% Linear SDDE definition
sysLDDE.F = @(t,p) p(1);    % present drift term
sysLDDE.Fd = @(t,p) p(2);   % delayed drift term
sysLDDE.f = @(t,p) p(3);    % constant drift term
sysLDDE.G = @(t,p) p(4);    % present diffusion term
sysLDDE.Gd = @(t,p) p(5);   % delayed diffusion term
sysLDDE.g = @(t,p) p(6);    % constant diffusion term


%% Find stationary correlation functions

na = 500; M = 10;
as = linspace(-5,5,na);
th = linspace(-tau,0);
phis = zeros(na,length(th));
eC = zeros(na,M); eF = zeros(1,na);
figure('Position',[400 200 1000 400])
for i = 1:na
    p_i = par; p_i(1) = as(i);
    phi_a = corr_scalar_asol(sysLDDE,tau,p_i);
    phis(i,:) = phi_a(th);
    JCi = pde_cov_jac(0,sysLDDE,p_i,tau,M);
    bCi = pde_cov_const(0,@(t)zeros(M,1),sysLDDE,p_i,tau,M);
    C_i = reshape(-JCi\reshape(bCi,[],1),M,M);
    eC(i,:) = eig(C_i);
    eF(i) = max(real(eig(JCi)));
end
eC(abs(eC)>10) = NaN; phis(abs(phis)>10) = NaN;
subplot(1,2,1); box on; hold on; grid on
surf(as,th,phis.','EdgeColor','none','FaceColor','b','FaceAlpha',0.3);
% surf(as,th,0*phis.','EdgeColor','none','FaceColor','k','FaceAlpha',0.5);
lightangle(45,45); light('Position',[-1 1 1]);
lighting gouraud; material shiny
C = orderedcolors("gem");
for i = 1:2:na
    plot3(as(i)+0*th,th,phis(i,:),'Color',C(1,:),'LineWidth',1)
end
plot3(as,zeros(1,na),phis(:,end),'k','LineWidth',2)
xlabel('$a$'); ylabel('$\vartheta$'); zlabel('$\bar{\phi}(\vartheta)$')
zlim([-3 5]); view(155,10)
annotation(gcf, 'textbox', [0.05, 0.9, 0.05, 0.05], ...
    'String', '(a)','LineStyle', 'none','Interpreter','latex');
fontsize(16,"points")
ax1 = subplot(2,2,2); box on; hold on; grid on; ax1.XTick = [-5 -2.5 0 2.5 5];
for i = 1:M
    plot(as,eC(:,i),'Color',(i-1)/M*[1 1 1],'LineWidth',2);
end
xline(-0.65625,':k');
ylim([-2.5 3]); xticklabels([]); ylabel('eig($\bar{C}(\theta,\vartheta)$)')
annotation(gcf, 'textbox', [0.47, 0.9, 0.05, 0.05], ...
    'String', '(b)','LineStyle', 'none','Interpreter','latex');
fontsize(16,"points")
ax2 = subplot(2,2,4); box on; hold on; grid on; ax2.XTick = [-5 -2.5 0 2.5 5];
pos = ax2.Position; pos(2) = pos(2) + 0.05; ax2.Position = pos;
mus = sqrt(as.^2-par(2)^2);
khi = ((par(4)+par(5))^2+2*(as+par(2))).*cosh(mus*par(7)/2) - ...
    (as+par(2))./mus.*((par(4)-par(5))^2+2*(as-par(2))).*sinh(mus*par(7)/2);
plot(as,khi,'r','LineWidth',2); hold on
plot(as,eF,'b','LineWidth',2);
yline(0,':k'); xline(-0.65625,':k'); hold off; axis([-5 5 -10 8]);
xlabel('$a$'); ylabel('$\chi$, $\rho(\mathbf{F})$')
annotation(gcf, 'textbox', [0.47, 0.5, 0.05, 0.05], ...
    'String', '(c)','LineStyle', 'none','Interpreter','latex');
 fontsize(16,"points")
exportgraphics(gcf,'plots/scalar_stat_ina1.png','Resolution',400)



