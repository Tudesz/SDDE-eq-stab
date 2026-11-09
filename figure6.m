clear; close all; clc;
%% Second moment stability of a linear, scalar SDDE
% dx(t) = (a*x(t)+b*x(t-tau) + c)*dt + 
%   (alpha*x(t)+beta*x(t-tau) + gamma)*dWt 
addpath(genpath('numerical methods'))

% COCO dependency
run("C:\Users\iklodi\Documents\MATLAB\coco_2025January28\startup.m");


%% System definition

% Default parameter vector
par = [0 -1 0 0 0 0 1]; % [a b c alpha beta gamma tau]
tau =@(p) p(7); d = 1;

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

% Stability map options
data.p0 = par;          % default parameter vector
data.pi_X = 5;          % index of parameter beta
data.pi_Y = 7;          % index of parameter tau     
data.num.print(2) = 0;  % no logging during generation
grid_res = [20 20];     % grid resolution
tau_lims = [1e-5 pi];   % parameter limits of tau
beta_lims = [0 2];      % parameter limits of beta

% COCO compatible functions wrappers
fc = @(f) @(p,d,u) deal(d,f(u));


%% First moment stability

data_m = data;
data_m.sys = sys_ps_approx_lode(sysLDDE,tau,d,2*M); 
timer = tic; 
fprintf('\nStability map generation:\n');
rest_m = zd_trianglestab(@condfunc_mean_stab, [beta_lims tau_lims],...
    2*grid_res, data_m, @f_stopcond);
rest_m.calc_time = toc(timer);
fprintf('   calculation time %0.3f s\n',rest_m.calc_time);
% figure();  plot_stab_lambda(rest_m,2);


%% Second moment stability DDE form

data_cD = data;
data_cD.sys = sysLDDE; 
data_cD.tau = tau; data_cD.M = M; 
timer = tic; 
fprintf('\nStability map generation:\n');
rest_cD = zd_trianglestab(@condfunc_cov_stab_del, [beta_lims tau_lims],...
    grid_res, data_cD, @f_stopcond);
rest_cD.calc_time = toc(timer);
fprintf('   calculation time %0.3f s\n',rest_cD.calc_time);
% figure();  plot_stab_lambda(rest_cD,2);


%% Save stability map data
if ~exist('data','dir')
    mkdir('data')
end
save('data/scalar_stab_map_atau.mat','rest_cD','rest_m')


%% "Analytic" stability boundary for the lambda = 0 case

% COCO compatible function wrapper (u = p0)
fch_eq = fc(@(u) sdde_char_sn(sysLDDE,u(7),u));
pnames = {'a' 'b' 'c' 'alpha' 'beta' 'gamma' 'tau'};
p0 = par; p0(1:2) = [0 -1];

% Create the COCO compatible continuation problem
prob = coco_prob();
prob = coco_add_func(prob,'ch_eq',fch_eq,[],'zero','u0',p0);
uidx = coco_get_func_data(prob, 'ch_eq', 'uidx');
prob = coco_add_pars(prob,'pars',uidx,pnames);

% Run COCO with b fixed
bd1 = coco(prob,'lsdde_at_psn1',[],1,{'tau','beta'},{tau_lims, beta_lims});


%% Combined stability map
load('data/scalar_stab_map_atau.mat','rest_cD','rest_m')

thm1.lspec = {'g-', 'LineWidth', 2};
figure('Position',[400 200 500 400]); 
plot_stab_lambda(rest_cD,2); hold on
contour = findobj(gca, 'Type', 'contour');
colormap('gray')
contour(2).LevelList = -4:0.5:2;
delete(contour(1))
c = findobj(gcf, 'Type', 'ColorBar');
c.Label.String = '$\rho(\mathbf{F})$';
c.TickLabelInterpreter = 'Latex';
plot_stab_lim(rest_m,1,{'k','LineWidth',2'});
coco_plot_bd(thm1,'lsdde_at_psn1','beta','tau');
bs = linspace(0,sqrt(2)); ts1 = 1-sqrt(bs.^2-bs.^4/4); % Eq.(67) from (Haskovec 2016)
ts2 = -bs.^2/2 + 1/4*sqrt(4*bs.^4+2*(2-bs.^2).^2); % Eq.(68) from (Haskovec 2016)
plot(bs,ts2,'g--',bs,ts1,'r--','LineWidth',2);
xlabel('$\beta$'); ylabel('$\tau$')
axis([0 2 0 2]); fontsize(16,"points")
exportgraphics(gcf,'plots/scalar_map_atau.png','Resolution',400)

