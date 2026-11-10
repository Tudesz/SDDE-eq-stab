clear; close all; clc;
%% Comparison of second moment stability maps of a linear, scalar SDDE
% dx(t) = (a*x(t)+b*x(t-tau) + c)*dt + 
%   (alpha*x(t)+beta*x(t-tau) + gamma)*dWt 
addpath(genpath('numerical methods'))


%% Analytic stability bounds

% Default parameter vector
par = [0 0 0 -1.5 0.0 0 1]; % [a b c alpha beta gamma tau]
tau = par(7);       % time delay
a_lims = [-4 1.5];  % parameter limits of a
b_lims = [-4 4];    % parameter limits of b

% Linear SDDE definition
sysLDDE.F = @(t,p) p(1);    % present drift term
sysLDDE.Fd = @(t,p) p(2);   % delayed drift term
sysLDDE.f = @(t,p) p(3);    % constant drift term
sysLDDE.G = @(t,p) p(4);    % present diffusion term
sysLDDE.Gd = @(t,p) p(5);   % delayed diffusion term
sysLDDE.g = @(t,p) p(6);    % constant diffusion term

% First-moment stability
n = 100;  % contour plot resolution
M = 10;   % discretization resolution 
a = linspace(a_lims(1),a_lims(2),n);
b = linspace(b_lims(1),b_lims(2),n);
ab_hopf = cosh(sqrt(a.^2-b.^2.')*tau/2)-...
    (a-b.')./sqrt(a.^2-b.^2.').*sinh(sqrt(a.^2-b.^2.')*tau/2);
ab_hopf(isnan(ab_hopf)) = 0;

% Second-moment stability
Ef = zeros(length(b),length(a));  % discretization-based
Ceq = zeros(length(b),length(a)); % analytic
for i = 1:size(Ceq,1)
    for j = 1:size(Ceq,2)
        p_temp = par; p_temp(1:2) = [a(j) b(i)];
        Ceq(i,j) = scalar_char_sn(sysLDDE,tau,p_temp);
        JCi = pde_cov_jac(0,sysLDDE,p_temp,tau,M);
        Ef(i,j) = max(real(eig(JCi)));
    end
end

% Erroneous stability predictions (Mackey 1995)
om = linspace(0,pi); 
hx = om.*cos(om*tau)./sin(om*tau);
hy = -om./sin(om*tau);
snx = linspace(a_lims(1),max(hx));
sny = -snx;


%% Plot results

figure('Position',[400 200 500 400]); 
contourf(a,b,Ef,[-inf -5:1:5],'LineStyle','none'); hold on;
c = colorbar; c.Label.Interpreter = 'Latex';
colormap('gray')
c.Label.String = '$\rho(\mathbf{F})$';
c.TickLabelInterpreter = 'Latex';
contour(a,b,Ceq,[0 0],'LineColor','g','LineWidth',2); 
plot(hx-1/2*par(4)^2,hy,'g--',snx-1/2*par(4)^2,sny,'g--','LineWidth',2)
plot(hx,hy,'k',snx,sny,'k','LineWidth',2)
axis([a_lims b_lims]); xlabel('$a\tau$'); ylabel('$b\tau$'); 
% annotation(gcf, 'textbox', [0.0, 0.9, 0.05, 0.05], ...
%     'String', '(a)','LineStyle', 'none','Interpreter','latex')
box on; fontsize(16,"points");
text(-2.95,0.95,'A','Color','w');
scatter(-3,1,'wx','LineWidth',2); 
text(-2.95,2.05,'B','Color','w');
scatter(-3,2.1,'wx','LineWidth',2); 
text(-0.35,0.-0.55,'C','Color','k');
scatter(-0.4,-0.9,'kx','LineWidth',2); 
text(1.05,0.-1.75,'D','Color','k');
scatter(1,-1.5,'kx','LineWidth',2); 
fontsize(16,'points')
% exportgraphics(gcf,'plots/scalar_stab_comp.png','Resolution',400)

