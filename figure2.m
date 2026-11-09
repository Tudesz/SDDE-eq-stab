clear; close all; clc;
%% First-moment stability of a linear, scalar SDDE
% dx(t) = (a*x(t)+b*x(t-tau) + c)*dt + 
%   (alpha*x(t)+beta*x(t-tau) + gamma)*dWt 
addpath(genpath('numerical methods'))

% RUN "scalar_moment_stab.m" BEFORE THIS SCIPT!

% COCO dependency
run("C:\Users\iklodi\Documents\MATLAB\coco_2025January28\startup.m");


%% System definition

% Default parameter vector
par = [0 0 0 -1.5 0.5 1 1]; % [a b c alpha beta gamma tau]
tau = par(7); d = 1;

% stochastic system definition
sysDDE.f = @(x,xd,t,p) p(1)*x(1,:) + p(2)*xd(1,:,1) + p(3); % drift term
sysDDE.g = @(x,xd,t,p) p(4)*x(1,:) + p(5)*xd(1,:,1) + p(6); % diffustion term

% Linear SDDE definition
sysLDDE.F = @(t,p) p(1);    % present drift term
sysLDDE.Fd = @(t,p) p(2);   % delayed drift term
sysLDDE.f = @(t,p) p(3);    % constant drift term
sysLDDE.G = @(t,p) p(4);    % present diffusion term
sysLDDE.Gd = @(t,p) p(5);   % delayed diffusion term
sysLDDE.g = @(t,p) p(6);    % constant diffusion term

% Discretization
M = 10; % Chebyshev mesh dimensions
a_lims = [-15 15];  % parameter limits of a
b_lims = [-15 15];  % parameter limits of b

% COCO compatible functions wrappers
fc = @(f) @(p,d,u) deal(d,f(u));


%% "Analytic" stability boundary for the lambda = 0 case

% COCO compatible function wrapper (u = p0)
fch_eq = fc(@(u) sdde_char_sn(sysLDDE,tau,u));
pnames = {'a' 'b' 'c' 'alpha' 'beta' 'gamma' 'tau'};
p0 = par; p0(1:2) = [-1 0];

% Create the COCO compatible continuation problem
prob = coco_prob();
prob = coco_add_func(prob,'ch_eq',fch_eq,[],'zero','u0',p0);
uidx = coco_get_func_data(prob, 'ch_eq', 'uidx');
prob = coco_add_pars(prob,'pars',uidx,pnames);

% Run COCO with b fixed
bd1 = coco(prob,'lsdde_psn1',[],1,{'a','b'},a_lims);


%% Closed form stability boundaries

tau = 1; data0.p0 = [0 0 2 -1.5 0.5 1 tau];
data0.pi_X = 1; data0.pi_Y = 2; data0.num.print(2) = 0; 
data0.sys = sys_ps_approx_lode(sysLDDE,tau,d,20); 
rest_m0 = zd_trianglestab(@condfunc_mean_stab, [a_lims b_lims],...
    [40 40], data0, @f_stopcond);
a = linspace(a_lims(1),a_lims(2),1000);
b = linspace(b_lims(1),b_lims(2),1000);
ab_hopf = cosh(sqrt(a.^2-b.^2.')*tau/2)-...
    (a-b.')./sqrt(a.^2-b.^2.').*sinh(sqrt(a.^2-b.^2.')*tau/2);
ab_hopf(isnan(ab_hopf)) = 0;
ab_fcond = cosh(sqrt(a.^2-b.^2.')*tau/2)-...
    (a+b.')./sqrt(a.^2-b.^2.').*sinh(sqrt(a.^2-b.^2.')*tau/2);
ab_fcond(isnan(ab_fcond)) = 0;
Ceq = zeros(length(b),length(a));
for i = 1:size(Ceq,1)
    for j = 1:size(Ceq,2)
        p_temp = par; p_temp(1:2) = [a(j) b(i)];
        Ceq(i,j) = scalar_char_sn(sysLDDE,tau,p_temp);
    end
end


%%
load('data/scalar_stab_map.mat','rest_m')
figure('Position',[400 200 600 400]); ax1 = gca;
plot_stab_lambda(rest_m,1); hold on
colormap([0.8 0.8 0.8; 1.0 1.0 1.0])
h = findobj(ax1, 'Type', 'Contour');
set(h,'EdgeColor','none','HandleVisibility','off')
patch(NaN, NaN, [0.8 0.2 0.2],'FaceColor','black','EdgeColor','none',...
    'FaceAlpha',0.2);
patch(NaN, NaN, [0.8 0.2 0.2],'FaceColor','blue','EdgeColor','none',...
    'FaceAlpha',20/255);
[~, h1] = contourf(a,b,ab_fcond,[0 inf],'EdgeColor','none',...
    'HandleVisibility','off');
% [~, h2] = contourf(a,b,Ceq,[0 inf],'EdgeColor','none',...
%     'HandleVisibility','off');
plot(a,-b,'k--','LineWidth',1.5);
plot(a,b,'k:','LineWidth',1.5,'HandleVisibility','off');
contour(a,b,ab_hopf,[0 0],'k','LineWidth',1.5); 
contour(a,b,ab_fcond,[0 0],'b','LineWidth',1.5); 
% contour(a,b,Ceq,[0 0],'r','LineWidth',1.5); 
axis([-15 15 -15 15]);  axis equal; hold off
xlabel('$a\tau$'); ylabel('$b\tau$');
legend('1$^{\mathrm{st}}$ m. stab.', '$\eta_{\tau/2}(a,-b)>0$','$a+b = 0$',...
    '$\eta_{\tau/2}(a,b) = 0$','$\eta_{\tau/2}(a,-b) = 0$','Location','northeastoutside'); 
% patch(NaN, NaN, [0.8 0.2 0.2],'FaceColor','red','EdgeColor','none',...
%     'FaceAlpha',20/255,'DisplayName','$\chi<0$'); 
fontsize(16,"points"); drawnow

face1 = h1.FacePrims(1);
face1.ColorType = 'truecoloralpha';
face1.ColorData = uint8([0; 0; 255; 20]); 
% face2 = h2.FacePrims(1);
% face2.ColorType = 'truecoloralpha';
% face2.ColorData = uint8([255; 0; 0; 20]);
exportgraphics(gcf,'plots/scalar_map_anal.png','Resolution',400)

