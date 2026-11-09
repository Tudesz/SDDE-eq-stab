clear; close all; clc;
%% Evaluate the second moment stability of a delayed oscillator
% dx(t) = (A*x(t) + B*x(t-tau))*dt + (alpha*x(t)+beta*x(t-tau)+gamma)*dW(t)
% fodor2023: A = [0 1; k 0]; B = [0 0; -p -d];
%   alpha = [0 0; 0 0]; beta = -sigma*B; gamma = [0; sigma]
% compare to Buckwar2013 for stability analysis
addpath(genpath('numerical methods'))

% COCO dependency
run("C:\Users\iklodi\Documents\MATLAB\coco_2025January28\startup.m");

%% System definition

% System matrices (fodor2023)
A = @(p) [0 1; p(4) 0];
B = @(p) [0 0; -p(1) -p(2)];
alpha = @(p) 0.0*[1 0; 0 -1];
beta = @(p) [0 0; p(3)*p(1) p(3)*p(2)];
gamma = @(p) [0; p(3)];

% Default parameter vector
par = [6 3 0.1 5 0.3]; % [p d sigma k tau]
tau = par(5); d = 2;

% stochastic system definition
sysDDE.f = @(x,xd,t,p) A(p)*x + B(p)*squeeze(xd(:,:,1));  % drift term
sysDDE.g = @(x,xd,t,p) alpha(p)*x + beta(p)*squeeze(xd(:,:,1)) + gamma(p); % diffustion term

% Linear SDDE definition
sysLDDE.F = @(t,p) A(p);        % present drift term
sysLDDE.Fd = @(t,p) B(p);       % delayed drift term
sysLDDE.f = @(t,p) [0; 0];      % constant drift term
sysLDDE.G = @(t,p) alpha(p);    % present diffusion term
sysLDDE.Gd = @(t,p) beta(p);    % delayed diffusion term
sysLDDE.g = @(t,p) gamma(p);    % constant diffusion term

% Discretization
M = 10; % Chebyshev mesh dimensions

% Stability map options
data.p0 = par;          % default parameter vector
data.pi_X = 1;          % index of parameter p
data.pi_Y = 2;          % index of parameter d     
data.num.print(2) = 0;  % no logging during generation
grid_res = [20 20];     % grid resolution
p_lims = [4 10];        % parameter limits of p
d_lims = [1 6];         % parameter limits of d


% COCO compatible functions wrappers
fc = @(f) @(p,d,u) deal(d,f(u));


%% First moment stability
data_m = data;
data_m.sys = sys_ps_approx_lode(sysLDDE,tau,d,2*M(end)); 
timer = tic; 
fprintf('\nStability map generation:\n');
rest_m = zd_trianglestab(@condfunc_mean_stab, [p_lims d_lims],...
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
rest_cD = zd_trianglestab(@condfunc_cov_stab_del, [p_lims d_lims],...
    grid_res, data_cD, @f_stopcond);
rest_cD.calc_time = toc(timer);
fprintf('   calculation time %0.3f s\n',rest_cD.calc_time);
% figure();  plot_stab_lambda(rest_cD{i},2);

% Save stability map data
if ~exist('data','dir')
    mkdir('data')
end
save('data/osc_stab_map.mat','rest_cD','rest_m')


%% Track the lambda = 0 stability boundary using continuation

% COCO compatible function wrapper (u = p0)
fch_eq = fc(@(u) sdde_char_sn(sysLDDE,tau,u));
pnames = {'p' 'd' 'sigma' 'k' 'tau'};

% Create the COCO compatible continuation problem
prob = coco_prob();
prob = coco_add_func(prob,'ch_eq',fch_eq,[],'zero','u0',par);
uidx = coco_get_func_data(prob, 'ch_eq', 'uidx');
prob = coco_add_pars(prob,'pars',uidx,pnames);

% Run COCO with both p and d released
prob = coco_set(prob,'cont','ItMX',[0 100]);
bd1 = coco(prob,'sdosc_psn1',[],1,{'p','d'},p_lims);

%% Form a surrogate DDE to find second moment stability (Buckwar 2013)

% Second moment characteristic equation
[Q, R, iu, il] = vecmat_symm(2);
P = vecmat_transp(d);
S = @(p) kron(eye(d),A(p)) + kron(A(p),eye(d)) ...
    + kron(alpha(p),alpha(p)); 
% T1 = @(p) (eye(d^2) + P) * (kron(B(p),eye(d)) + kron(beta(p),alpha(p)));
T1 = @(p) kron(eye(d),B(p)) + kron(B(p),eye(d)) + ...
    + kron(beta(p),alpha(p)) + kron(alpha(p),beta(p));
T2 = @(p) kron(beta(p),beta(p));
% char_eq = @(lam,p,tau) det(lam*eye(4^2) - 1/2*S(p) - 1/2*T1(p)*exp(-lam*tau)...
%     -1/2*T2(p)*exp(-2*lam*tau)); 

% Surrogate DDE for stability analysis 
tau2 = [tau, 2*tau];                % extended delay vector
sys2LDDE.F = @(t,p) 1/2*S(p);       % present drift term
sys2LDDE.Fd{1} = @(t,p) 1/2*T1(p);  % delayed drift term
sys2LDDE.Fd{2} = @(t,p) 1/2*T2(p);  % delayed drift term
sys2LDDE.f = @(t,p) zeros(d^2,1);   % constant drift term
sys2LDDE.G = @(t,p) zeros(d^2);     % present diffusion term
sys2LDDE.Gd{1} = @(t,p) zeros(d^2); % delayed diffusion term
sys2LDDE.Gd{2} = @(t,p) zeros(d^2); % delayed diffusion term
sys2LDDE.g = @(t,p) zeros(d^2,1);   % constant diffusion term

% "First" moment stability of the surragate system
data_m = data;
data_m.sys = sys_ps_approx_lode(sys2LDDE,tau2,d^2,M); 
timer = tic; 
fprintf('\nStability map generation:\n');
rest_m2 = zd_trianglestab(@condfunc_mean_stab, [p_lims d_lims],...
    2*grid_res, data_m, @f_stopcond);
rest_m2.calc_time = toc(timer);
fprintf('   calculation time %0.3f s\n',rest_m2.calc_time);
% figure(); plot_stab_lambda(rest_m2,2)

% save data
save('data/osc_stab_map_comp.mat','rest_m2')


%% Combined stability map
load('data/osc_stab_map.mat','rest_m','rest_cD')
load('data/osc_stab_map_comp.mat','rest_m2')

thm1.lspec = {'g-', 'LineWidth', 1.5};
thm2.lspec = {'--', 'LineWidth', 1.5};
% figure('Position',[400 400 400 400]);
figure('Position',[400 400 600 500]);
ax1 = axes('Position',[0.1138 0.1482 0.77 0.7768]); box on;
plot_stab_lambda(rest_cD,2); hold on
contour = findobj(gca, 'Type', 'contour');
colormap('gray')
contour(2).LevelList = -3:0.5:3;
delete(contour(1))
c = findobj(gcf, 'Type', 'ColorBar');
c.Label.String = '$\rho(\mathbf{F})$'; hold on
c.TickLabelInterpreter = 'Latex';
set(get(gca, 'Children'), 'HandleVisibility', 'off'); 
plot_stab_lim(rest_m,1,{'k','LineWidth',1.5'});
plot_stab_lim(rest_m2,1,{'g--','LineWidth',1.5'});
coco_plot_bd(thm1,'sdosc_psn1','p','d');
scatter(6.5,3.5,'wx','LineWidth',2); 
text(6.55,3.65,'A','Color','w','FontSize',16);
xlabel('$p$'); ylabel('$d$'); axis([4 10 1 6]);  hold off
% legend('first moment','second moment','Location','southeast')
ax1.Units = 'normalized'; pos = ax1.Position; 
annotation(gcf, 'textbox', [0, pos(2)+pos(4), 0.04, 0.04], ...
    'String', '(d)','LineStyle', 'none','Interpreter','latex',...
    'HorizontalAlignment','left');
fontsize(20,"points")
exportgraphics(gcf,'plots/osc_map_comp.png','Resolution',400)