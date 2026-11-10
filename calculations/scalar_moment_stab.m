clear; close all; clc;
%% First and second-moment stability maps of a linear, scalar SDDE
% dx(t) = (a*x(t)+b*x(t-tau) + c)*dt + 
%   (alpha*x(t)+beta*x(t-tau) + gamma)*dWt 
addpath(genpath('../numerical methods'))


%% System definition

% Default parameter vector
par = [0 0 0 1 0 1 1];  % [a b c alpha beta gamma tau]
tau = par(7);           % time delay
a_lims = [-15 15];      % parameter limits of a
b_lims = [-15 15];      % parameter limits of b

% Linear SDDE definition
sysLDDE.F = @(t,p) p(1);    % present drift term
sysLDDE.Fd = @(t,p) p(2);   % delayed drift term
sysLDDE.f = @(t,p) p(3);    % constant drift term
sysLDDE.G = @(t,p) p(4);    % present diffusion term
sysLDDE.Gd = @(t,p) p(5);   % delayed diffusion term
sysLDDE.g = @(t,p) p(6);    % constant diffusion term

% Stability map options
M = 10;                 % discretization dimension
data.p0 = par;          % default parameter vector
data.pi_X = 1;          % index of parameter a
data.pi_Y = 2;          % index of parameter b     
data.num.print(2) = 0;  % no logging during generation
grid_res = [5 5];       % grid resolution


%% First-moment stability
data_m = data;
data_m.sys = sys_ps_approx_lode(sysLDDE,tau,1,2*M); 
timer = tic; 
fprintf('\nStability map generation:\n');
rest_m = zd_trianglestab(@condfunc_mean_stab, [a_lims b_lims],...
    2*grid_res, data_m, @f_stopcond);
rest_m.calc_time = toc(timer);
fprintf('   calculation time %0.3f s\n',rest_m.calc_time);
figure();  plot_stab_lambda(rest_m,2);



%% Second- moment stability
data_cD = data;
data_cD.sys = sysLDDE; 
data_cD.tau = tau; data_cD.M = M; 
timer = tic; 
fprintf('\nStability map generation:\n');
rest_cD = zd_trianglestab(@condfunc_cov_stab_del, [a_lims b_lims],...
    grid_res, data_cD, @f_stopcond);
rest_cD.calc_time = toc(timer);
fprintf('   calculation time %0.3f s\n',rest_cD.calc_time);
figure();  plot_stab_lambda(rest_cD,2);


%% Export data
save('../data/scalar_first_stab_map.mat','rest_m','rest_cD')

