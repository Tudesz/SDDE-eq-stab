clear; close all; clc;
%% Evaluate the second moment stability of a stochastic turning operation
% dx(t) = (a*x(t) + b*x(t-tau))*dt + (alpha*x(t)+beta*x(t-tau)+gamma)*dW(t)
% a = [0 1; -1-w -2*zeta]; b = [0 0; w 0]; alpha = [0 0; -sigma 0]; 
% beta = [0 0; sigma 0]; gamma = [0; sigma]; tau = 2*pi/omega
addpath(genpath('numerical methods'))

% COCO dependency
run("C:\Users\iklodi\Documents\MATLAB\coco_2025January28\startup.m");


%% System definition

% System matrices
A = @(p) [0 1; -1-p(2) -2*p(1)];
B = @(p) [0 0; p(2) 0];
alpha = @(p) [0 0; -p(4)*p(2) 0];
beta = @(p) [0 0; p(4)*p(2) 0];
gamma = @(p) [0; p(4)*p(2)];

% Default parameter vector
par = [0.05 0.1 2 1.0]; % [zeta w omega sigma]
tau = 2*pi/par(3); % time delay

% Linear SDDE definition
sysLDDE.F = @(t,p) A(p);            % present drift term
sysLDDE.Fd{1} = @(t,p) B(p);        % delayed drift term
sysLDDE.f = @(t,p) [0; 0];          % constant drift term
sysLDDE.G = @(t,p) alpha(p);        % present diffusion term
sysLDDE.Gd{1} = @(t,p) beta(p);     % delayed diffusion term
sysLDDE.g = @(t,p) gamma(p);        % constant diffusion term

% Stability map options
data.p0 = par;          % default parameter vector
data.pi_X = 3;          % index of parameter om
data.pi_Y = 2;          % index of parameter w    
data.num.print(2) = 0;  % no logging during generation
grid_res = [20 20];     % grid resolution
om_lims1 = [0.1 0.4];   % parameter limits of om
om_lims2 = [0.4 2.0];   % parameter limits of om
w_lims = [0 0.8];       % parameter limits of w


%% Second moment stability from stationary correlation positive definiteness

data_cO = data; data_cO.tau = @(p) 2*pi/p(3);
data_cO.sys = sysLDDE;
timer = tic; 
fprintf('\nStability map generation:\n');
rest_cO1 = zd_trianglestab(@condfunc_stat_stab,[om_lims1 w_lims],...
    grid_res, data_cO, @f_stopcond);
 rest_cO1.calc_time = toc(timer);
fprintf('   calculation time %0.3f s\n',rest_cO1.calc_time);
% figure(); plot_stab_lambda(rest_cO1,2);
% figure(); plot_stab_test(rest_cO1);
fprintf('\nStability map generation:\n');
rest_cO2 = zd_trianglestab(@condfunc_stat_stab,[om_lims2 w_lims],...
grid_res, data_cO, @f_stopcond);
 rest_cO2.calc_time = toc(timer);
fprintf('   calculation time %0.3f s\n',rest_cO2.calc_time);
% figure(); plot_stab_lambda(rest_cO2,2);
% figure(); plot_stab_test(rest_cO2);


%% Save stability map data
if ~exist('data','dir')
    mkdir('data')
end
save('data/sturn_stab_map_a.mat','rest_cO1','rest_cO2')

%% Plot results
load('data/sturn_stab_map_a.mat','rest_cO1','rest_cO2')
thm1.lspec = {'r-', 'LineWidth', 1.5};
figure('Position',[400 400 800 350]); box on;
plot_stab_lambda(rest_cO1,1);  hold on;
plot_stab_lambda(rest_cO2,1);  hold on;
coco_plot_bd(thm1,'sturn_psn1','omega','w'); hold off
xlim([0.1 2.0])

% figure(); 
% subplot(2,1,1); plot_stab_test(rest_cO1,1);
% subplot(2,1,2); plot_stab_test(rest_cO1,2);


