clear; close all; clc;
%% Evaluate the stationary variance of a stochastic turning operation
% dx(t) = (a*x(t) + b*x(t-tau))*dt + (alpha*x(t)+beta*x(t-tau)+gamma)*dW(t)
% a = [0 1; -1-w -2*zeta]; b = [0 0; w 0]; alpha = [0 0; -sigma 0]; 
% beta = [0 0; sigma 0]; gamma = [0; sigma]; tau = 2*pi/omega
addpath(genpath('numerical methods'))


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

% Discretization
M = 40; % Chebyshev mesh dimensions

% Stability map options
data.p0 = par;          % default parameter vector
data.pi_X = 3;          % index of parameter om
data.pi_Y = 2;          % index of parameter w    
data.num.print(2) = 0;  % no logging during generation
grid_res = [20 20];     % grid resolution
om_lims = [0.1 2];      % parameter limits of om
w_lims = [0 0.8];       % parameter limits of w

% COCO compatible functions wrappers
fc = @(f) @(p,d,u) deal(d,f(u));
f_tau = fc(@(u) u(5)-2*pi/u(3));
fch_eq = fc(@(u) sdde_char_sn(sysLDDE,u(5),u(1:4)));
pnames = {'zeta', 'w', 'omega', 'sigma', 'tau'};


%% First-moment stability boundary

% Find branches of Hopf bifurcations
nl = 15; % number of lobes
om_c = linspace(1+1e-3,2,1000); % chatter frequency
zeta = par(1);
wD = ((om_c.^2-1).^2+4*zeta^2*om_c.^2)./(2*(om_c.^2-1));
omD = zeros(nl,length(om_c));
for k = 1:nl
    omD(k,:) = pi*om_c./(atan((1-om_c.^2)./(2*zeta*om_c))+k*pi);
end

% Find the lowest lobe / actual stability boundary
omA = linspace(om_lims(1),om_lims(2),2000);
wA = inf*ones(size(omA));
for i = 1:length(omA)
    for j = find(omD(:,1)<omA(i) & omA(i)<omD(:,end)).'
        w_temp = interp1(omD(j,:),wD,omA(i));
        if w_temp < wA(i)
            wA(i) = w_temp;
        end
    end
end
 % plot(omA,wA);


 %% Second moment stability boundary for the lambda = 0 case

% Initial guess for lambda ~ 0
p1 = par; p1(2:3) = [0.1 1.3]; tau1 = 2*pi/p1(3); u0 = [p1 tau1];

% Create the COCO compatible continuation problem
prob = coco_prob();
prob = coco_add_func(prob,'tau_eq',f_tau,[],'zero','u0',u0);
uidx = coco_get_func_data(prob, 'tau_eq', 'uidx');
prob = coco_add_func(prob,'ch_eq',fch_eq,[],'zero','uidx',uidx);
prob = coco_add_pars(prob,'pars',uidx,pnames);

% Run COCO
prob = coco_set(prob,'cont','norm',inf);
% prob = coco_set(prob,'cont','h_max',0.1);
prob = coco_set(prob,'cont','ItMX',1000);
bd1 = coco(prob,'sturn_psn1',[],1,{'w','omega','tau'},{w_lims,om_lims});


%% Stationary variance heatmap
[om2,iom] = unique(coco_bd_col(bd1,'omega')); 
w2 = coco_bd_col(bd1,'w'); w2 = w2(iom);
na = 1000;
oms = linspace(om_lims(1),om_lims(2),na);
ws = linspace(w_lims(1),w_lims(2),na);
Vx = NaN(na,na);
fprintf('\nContourplot generation\n');
for i = 1:na
    taui = 2*pi/oms(i);
    for j = 1:na
        pij = par; pij(2:3) = [ws(j) oms(i)];
        if ws(j) < interp1(om2,w2,oms(i))
            phi = corr_stat_asol2(sysLDDE,taui,pij); 
            f0_ij = phi(0);
            Vx(j,i) = f0_ij(1,1);
        end
    end
    fprintf('  progress: %i/%i\n',i,na)
end

%% Plot stability map
thm1.lspec = {'r-', 'LineWidth', 1.5};
clrmp = orderedcolors("gem");
figure('Position',[400 400 600 400]); box on; hold on;
vxl = log10([0 logspace(-2,log10(20*par(4)^2),20)]);
% vxl = log10(linspace(0,5*par(4)^2,20));
contourf(oms,ws,log10(abs(Vx)),vxl,'EdgeColor','w','LineWidth',0.1,...
    'LineStyle','none','HandleVisibility','off','FaceAlpha',1.0);
c = colorbar; colormap("jet")
c.Label.Interpreter = 'Latex';
c.TickLabelInterpreter = 'Latex';
c.Label.String = 'log$_{10}$($\bar{\phi}_{11}(0)$)';
% plot(omA, wA, 'k', 'LineWidth', 1.5);
coco_plot_bd(thm1,'sturn_psn1','omega','w');
xlabel('$\tilde{\omega}$'); ylabel('$\tilde{w}$'); hold off
% legend('mean stab.', 'mean-square stab.'); 
annotation(gcf, 'textbox', [0.0, 0.92, 0.05, 0.05], ...
    'String', '(c)','LineStyle', 'none','Interpreter','latex');
title('$\zeta = 0.05$, $\sigma = 1$')
fontsize(16,"points"); axis([0.1 2 0 0.57])
exportgraphics(gcf,'plots/turning_map_var.png','Resolution',400)
