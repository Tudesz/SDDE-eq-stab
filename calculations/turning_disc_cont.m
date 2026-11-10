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

% Parameter Jacobians
A_dp = @(p) cat(3,[0 0; 0 -2], [0 0; -1 0], zeros(2), zeros(2));
B_dp = @(p) cat(3,zeros(2), [0 0; 1 0], zeros(2), zeros(2));
alpha_dp = @(p) cat(3,zeros(2), [0 0; -p(4) 0], zeros(2), [0 0; -p(2) 0]);
beta_dp = @(p) cat(3,zeros(2), [0 0; p(4) 0], zeros(2), [0 0; p(2) 0]);
gamma_dp = @(p) cat(3,zeros(2,1), [0; p(4)], zeros(2,1), [0; p(2)]);
tau_dp = @(p) cat(3, 0, 0, -2*pi/p(3)^2, 0);

% stochastic system definition
sysDDE.f = @(x,xd,t,p) A(p)*x + B(p)*squeeze(xd(:,:,1));  % drift term
sysDDE.g = @(x,xd,t,p) alpha(p)*x + beta(p)*squeeze(xd(:,:,1)) + gamma(p); % diffustion term

% Linear SDDE definition
sysLDDE.F = @(t,p) A(p);            % present drift term
sysLDDE.Fd{1} = @(t,p) B(p);        % delayed drift term
sysLDDE.f = @(t,p) [0; 0];          % constant drift term
sysLDDE.G = @(t,p) alpha(p);        % present diffusion term
sysLDDE.Gd{1} = @(t,p) beta(p);     % delayed diffusion term
sysLDDE.g = @(t,p) gamma(p);        % constant diffusion term


% Linear SDDE parameter Jacobians
sel_pi = @(x,i) x(:,:,i);
sysLDDE.F_p = @(t,p,i) sel_pi(A_dp(p),i);           % present drift term
sysLDDE.Fd_p{1} = @(t,p,i) sel_pi(B_dp(p),i);       % delayed drift term
sysLDDE.f_p = @(t,p,i) [0; 0];                      % constant drift term
sysLDDE.G_p = @(t,p,i) sel_pi(alpha_dp(p),i);       % present diffusion term
sysLDDE.Gd_p{1} = @(t,p,i) sel_pi(beta_dp(p),i);    % delayed diffusion term
sysLDDE.g_p = @(t,p,i) sel_pi(gamma_dp(p),i);       % constant diffusion term
sysLDDE.tau_p = @(t,p,i) sel_pi(tau_dp(p),i);       % time delay vector

% Discretization
M = [5 10 20 40]; % Chebyshev mesh dimensions

% Stability map options
data.p0 = par;          % default parameter vector
data.pi_X = 3;          % index of parameter om
data.pi_Y = 2;          % index of parameter w    
data.num.print(2) = 0;  % no logging during generation
grid_res = [20 20];     % grid resolution
om_lims = [0.1 2];      % parameter limits of om
w_lims = [0 0.8];       % parameter limits of w

% COCO compatible functions wrappers
pnames = {'zeta', 'w', 'omega', 'sigma'};
mesh.n = 2; mesh.tau = @(p) 2*pi/p(3);

%% Second moment stability lobes using collocation based continuation
% 2 variable correlation function BVP

for i = 1:length(M)

    % Initialization
    Mi = M(i); mesh.M = Mi;
    rname1 = sprintf('sturn_pp1_M%i',M(i));
    rname2 = sprintf('sturn_pp2_M%i',M(i));

    % Initial solution guess at a point where lambda ~ 0
    p2 = par; p2(2:3) = [0.11 0.6]; tau2 = 2*pi/p2(3);
    Jphi = pde_corr_jac(0,sysLDDE,p2,tau2,Mi);
    [V,D] = eig(Jphi); lamC = diag(D);
    ic = 1:length(lamC); ic(imag(lamC)~=0) = [];
    [lamc, ic0] = max(real(lamC(ic)));
    Vc = V(:,ic(ic0));

    % Create the COCO compatible continuation problem
    u0 = [lamc; Vc; p2.']; u_idx = 1:length(u0);
    par_idx = 2+length(Vc):length(u0);
    fcM_eq = @(p,d,u) corr_zero(p,d,u,sysLDDE,mesh);
    prob = coco_prob();
    % prob = coco_add_func(prob,'ch_eq',fcM_eq,[],'zero','u0',u0);
    prob = coco_add_func(prob,'ch_eq',fcM_eq,[],'zero','u0',u0,'f+df');
    prob = coco_add_pars(prob,'pars',[1 par_idx],[{'lambda'} pnames]);
    prob = coco_add_event(prob,'SNL','boundary','lambda','=',0);

    % Run COCO with d fixed
    bd3 = coco(prob,rname1,[],1,{'w','lambda'},[0 2]);

    % Run COCO again with a fixed lambda
    prob = coco_prob();
    lab = coco_bd_labs(rname1,'SNL');
    chart = coco_read_solution(rname1,lab(1),'chart');
    % prob = coco_add_func(prob,'ch_eq',fcM_eq,[],'zero',...
    %     'u0',chart.x(u_idx));
    prob = coco_add_func(prob,'ch_eq',fcM_eq,[],'zero',...
        'u0',chart.x(u_idx),'f+df');
    prob = coco_add_pars(prob,'pars',[1 par_idx],[{'lambda'} pnames]);
    prob = coco_set(prob,'cont','norm',inf);
    prob = coco_set(prob,'cont','ItMX',1000);
    rt = tic;
    bd4 = coco(prob,rname2,[],1,{'w','omega'},{w_lims, om_lims});
    toc(rt)
end


%% Second moment stability lobes using collocation based continuation
% 3 variable correlation function BVP
% 
% for i = 1:length(M)
%     % Initialization
%     Mi = M(i);
%     rname1 = sprintf('sturn_pc1_M%i',M(i));
%     rname2 = sprintf('sturn_pc2_M%i',M(i));
% 
%     % Initial solution guess at a point where lambda ~ 0
%     p2 = par; p2(2:3) = [0.11 0.6]; tau2 = 2*pi/p2(3);
%     [Qi, Ri, ~, ~] = vecmat_symm(Mi*2);
%     JC = pde_cov_jac(0,sysLDDE,p2,tau2,Mi);
%     % [V,D] = eig(JC);
%     [V,D] = eig(Qi*JC*Ri);
%      lamC = diag(D); [~,ic] = max(real(lamC)); 
%     Vc = V(:,ic); lamc = lamC(ic);
% 
%     % Create the COCO compatible continuation problem
%     u0 = [lamc; Vc; p2.']; u_idx = 1:length(u0);
%     par_idx = 2+length(Vc):length(u0);
%     mesh.M = Mi; mesh.Q = Qi; mesh.R = Ri;
%     fcM_eq = @(p,d,u) cov_red_zero(p,d,u,sysLDDE,mesh);
%     prob = coco_prob();
%     % prob = coco_add_func(prob,'ch_eq',fcM_eq,[],'zero','u0',u0);
%     prob = coco_add_func(prob,'ch_eq',fcM_eq,[],'zero','u0',u0,'f+df');
%     prob = coco_add_pars(prob,'pars',[1 par_idx],[{'lambda'} pnames]);
%     prob = coco_add_event(prob,'SNL','boundary','lambda','=',0);
% 
%     % Run COCO with d fixed
%     bd3 = coco(prob,rname1,[],1,{'w','lambda'},w_lims);
% 
%     % Run COCO again with a fixed lambda
%     prob = coco_prob();
%     lab = coco_bd_labs(rname1,'SNL');
%     chart = coco_read_solution(rname1,lab(1),'chart');
%     % prob = coco_add_func(prob,'ch_eq',fcM_eq,[],'zero',...
%     %     'u0',chart.x(u_idx));
%     prob = coco_add_func(prob,'ch_eq',fcM_eq,[],'zero',...
%         'u0',chart.x(u_idx),'f+df');
%     prob = coco_add_pars(prob,'pars',[1 par_idx],[{'lambda'} pnames]);
%     prob = coco_set(prob,'cont','norm',inf);
%     prob = coco_set(prob,'cont','ItMX',1000);
%     bd4 = coco(prob,rname2,[],1,{'w','omega'},{w_lims, om_lims});
% end


%% Plot stability map

thm1.lspec = {'r-', 'LineWidth', 1.5};
thm2.lspec = {'--', 'LineWidth', 1.5};
thm3.lspec = {':', 'LineWidth', 1.5};
clrmp = orderedcolors("gem");
figure('Position',[400 400 800 350]); box on; hold on;
coco_plot_bd(thm1,'sturn_psn1','omega','w');
for Mi = M
    coco_plot_bd(thm3,sprintf('sturn_pp2_M%i',Mi),'omega','w');
    coco_plot_bd(thm2,sprintf('sturn_pc2_M%i',Mi),'omega','w');
end
xlabel('$\tilde{\omega}$'); ylabel('$\tilde{w}$'); hold off
legend('first moment','second moment'); axis([0.1 2.0 0.0 0.8])
% % title(sprintf('$\\zeta$ = %0.2f, $\\sigma$ = %0.1f',par(1),par(4)))
fontsize(16,"points")
% exportgraphics(gcf,'plots/turning_map.png','Resolution',400)

