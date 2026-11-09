clear; close all; clc;
%% Evaluate the second-moment stability of a delayed oscillator
% dx(t) = (A*x(t) + B*x(t-tau))*dt + (alpha*x(t)+beta*x(t-tau)+gamma)*dW(t)
% fodor2023: A = [0 1; k 0]; B = [0 0; -p -d];
%   alpha = [0 0; 0 0]; beta = -sigma*B; gamma = [0; sigma]
addpath(genpath('numerical methods'))

% COCO dependency
run("C:\Users\iklodi\Documents\MATLAB\coco_2025January28\startup.m");


%% System definition

% System matrices (fodor2023)
A = @(p) [0 1; p(4) 0];
B = @(p) [0 0; -p(1) -p(2)];
beta = @(p) [0 0; p(3)*p(1) p(3)*p(2)];
gamma = @(p) [0; p(3)];

% Default parameter vector
par = [5.2 3 0.1 5 0.3]; % [p d sigma k tau]
tau = par(5); d = 2;

% stochastic system definition
sysDDE.f = @(x,xd,t,p) A(p)*x + B(p)*squeeze(xd(:,:,1));  % drift term
sysDDE.g = @(x,xd,t,p) alpha(p)*x + beta(p)*squeeze(xd(:,:,1)) + gamma(p); % diffustion term

% Linear SDDE definition
sysLDDE.F = @(t,p) A(p);        % present drift term
sysLDDE.Fd = @(t,p) B(p);       % delayed drift term
sysLDDE.f = @(t,p) [0; 0];      % constant drift term
sysLDDE.G = @(t,p) zeros(2);    % present diffusion term
sysLDDE.Gd = @(t,p) beta(p);    % delayed diffusion term
sysLDDE.g = @(t,p) gamma(p);    % constant diffusion term

% Discretization
M = [5 10]; % Chebyshev mesh dimensions

% Stability map options
p_lims = [5 15];      % parameter limits of p
d_lims = [1.25 6.55];     % parameter limits of d
tau_lims = [0.25 0.65];  % parameter limits of d

% COCO compatible functions wrappers
fc = @(f) @(p,d,u) deal(d,f(u));


%% Track the lambda = 0 stability boundary using continuation

% COCO compatible function wrapper (u = [lambda; p0])
fch_eq = @(~,d,u) deal(d,sdde_char_sn(sysLDDE,u(5),u));
pnames = {'p' 'd' 'sigma' 'k' 'tau'};

% Create the COCO compatible continuation problem
prob = coco_prob();
prob = coco_add_func(prob,'ch_eq',fch_eq,[],'zero','u0',par);
uidx = coco_get_func_data(prob, 'ch_eq', 'uidx');
prob = coco_add_pars(prob,'pars',uidx,pnames);

% Run COCO with both p and d released
prob = coco_set(prob, 'cont', 'atlas', 'kd','PtMX',5000,...
    'R', 0.1, 'R_max', 10, 'R_min', 1e-5,'al_max',15);
bd1 = coco(prob,'sdosc_sn1_kd',[],2,...
    {'p','d','tau'},{p_lims,d_lims,tau_lims});


%% Plot stability map
       
% x  = coco_bd_col(bd1, {'p' 'd' 'tau'});
% plot3(x(1,:), x(2,:), x(3,:), 'LineStyle', 'none', 'LineWidth', 2, ...
%   'Color', 'black', 'Marker', '.', 'MarkerSize', 15, ...
%   'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'black')

% 3d COCO output
figure('Position',[400 200 600 500]);
atlas = coco_bd_read('sdosc_sn1_kd', 'atlas');
plot_atlas_kd(atlas.charts,3); hold on
h = findobj(gca, 'Type', 'Patch');
set(h,'FaceColor','r','HandleVisibility','off','FaceAlpha',0.7)

% Analytic first moment stability boundary surface
tau = linspace(tau_lims(1),tau_lims(2),100);
oms = linspace(0,2*pi,100);
pas = (oms.^2+par(4)).*cos(tau.'*oms);
das = ((oms.^2+par(4))./oms).*sin(tau.'*oms);
taus = repmat(tau,length(oms),1);
surf(pas,das,taus.','LineStyle','none','FaceColor','b',...
    'FaceAlpha',0.2); hold on;
psn = par(4)*ones(size(pas));
dsn = repmat(linspace(d_lims(1),d_lims(2),length(taus)),length(oms),1);
surf(psn, dsn, taus.','LineStyle','none','FaceColor','k',...
    'FaceAlpha',0.2);
for i = [1:10:length(tau) length(tau)]
    plot3(pas(i,:),das(i,:),taus(:,i),'b','LineWidth',0.5,...
        'HandleVisibility','off')
    plot3(pas(:,i),das(:,i),taus(i,:),'b','LineWidth',0.5,...
        'HandleVisibility','off')
end
scatter3(par(4),sqrt(2*par(4)),sqrt(2/par(4)),50,'bo','filled')
xlim([5 12]); ylim([0 7]); zlim(tau_lims)
xlabel('$p$'); ylabel('$d$'); zlabel('$\tau$')
lightangle(45,45); light('Position',[-1 1 1]); view(60,30)
lighting gouraud; material shiny; box on; grid on;
annotation(gcf, 'textbox', [0, 0.9, 0.05, 0.05], ...
    'String', '(e)','LineStyle', 'none','Interpreter','latex');
fontsize(20,"points");
exportgraphics(gcf,'plots/stoch_osc_3D.png','Resolution',400)