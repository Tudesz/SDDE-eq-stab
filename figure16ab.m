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


%% First-moment stability boundary for different damping levels
zeta = [0.05 0.1 0.15 0.2];

% figure(); hold on;
omA = cell(1,length(zeta));
wA = cell(1,length(zeta));
for l = 1:length(zeta)
    % Find branches of Hopf bifurcations
    nl = 15; % number of lobes
    om_c = linspace(1+1e-3,2,1000); % chatter frequency
    wD = ((om_c.^2-1).^2+4*zeta(l)^2*om_c.^2)./(2*(om_c.^2-1));
    omD = zeros(nl,length(om_c));
    for k = 1:nl
        omD(k,:) = pi*om_c./(atan((1-om_c.^2)./(2*zeta(l)*om_c))+k*pi);
    end
    
    % Find the lowest lobe / actual stability boundary
    omA{l} = linspace(om_lims(1),om_lims(2),2000);
    wA{l} = inf*ones(size(omA{l}));
    for i = 1:length(omA{l})
        for j = find(omD(:,1)<omA{l}(i) & omA{l}(i)<omD(:,end)).'
            w_temp = interp1(omD(j,:),wD,omA{l}(i));
            if w_temp < wA{l}(i)
                wA{l}(i) = w_temp;
            end
        end
    end
     % plot(omA{l},wA{l}); drawnow
end


%% Second-moment stability boundary for different noise levels
sigma = [1 1.5 2 2.5];
rnames1 = cell(1,length(sigma));

for i = 1:length(sigma)
    % Initial guess for lambda ~ 0
    p1 = par; p1(2:3) = [0.1 1.3]; tau1 = 2*pi/p1(3); 
    p1(4) = sigma(i); u0 = [p1 tau1];

    % Create the COCO compatible continuation problem
    prob = coco_prob();
    prob = coco_add_func(prob,'tau_eq',f_tau,[],'zero','u0',u0);
    uidx = coco_get_func_data(prob, 'tau_eq', 'uidx');
    prob = coco_add_func(prob,'ch_eq',fch_eq,[],'zero','uidx',uidx);
    prob = coco_add_pars(prob,'pars',uidx,pnames);
    
    % Run COCO
    % prob = coco_set(prob,'cont','norm',inf);
    prob = coco_set(prob,'cont','h_max',0.25*sigma(i)/sigma(end));
    prob = coco_set(prob,'cont','ItMX',1000);
    rnames1{i} = sprintf('sturn_sigma%i_psn1',i);
    bd = coco(prob,rnames1{i},[],1,{'w','omega','tau'},{w_lims,om_lims});
end


%% Second-moment stability boundary for different damping levels
rnames2 = cell(1,length(zeta));

for i = 1:length(zeta)
    % Initial guess for lambda ~ 0
    p1 = par; p1(2:3) = [0.1*zeta(i)/par(1) 1.3]; tau1 = 2*pi/p1(3); 
    p1(1) = zeta(i); u0 = [p1 tau1];

    % Create the COCO compatible continuation problem
    prob = coco_prob();
    prob = coco_add_func(prob,'tau_eq',f_tau,[],'zero','u0',u0);
    uidx = coco_get_func_data(prob, 'tau_eq', 'uidx');
    prob = coco_add_func(prob,'ch_eq',fch_eq,[],'zero','uidx',uidx);
    prob = coco_add_pars(prob,'pars',uidx,pnames);
    
    % Run COCO
    % prob = coco_set(prob,'cont','norm',inf);
    if i > 1
        prob = coco_set(prob,'cont','h_min',0.2*zeta(i)/zeta(end));
    else
        prob = coco_set(prob,'cont','h_max',0.1);
    end
    prob = coco_set(prob,'cont','ItMX',1000);
    rnames2{i} = sprintf('sturn_zeta%i_psn1',i);
    bd = coco(prob,rnames2{i},[],1,{'w','omega','tau'},{w_lims,om_lims});
end


%% Plot stability maps
sigma = [1 1.5 2 2.5];
rnames1 = sprintfc('sturn_sigma%i_psn1',1:4);

clrmp = orderedcolors("gem");
figure('Position',[400 400 600 400]); box on; hold on;
a = area(omA{1}, wA{1}); a.HandleVisibility = 'off';
a.FaceColor = [0.9 0.9 0.9]; a.EdgeColor = 'none';
ax = gca; ax.Layer = 'top'; 
plot(omA{1}, wA{1}, 'k', 'LineWidth', 1.5);
for i = 1:length(sigma)
    thmi.lspec = {'-', 'LineWidth', 1.5,'Color',clrmp(i,:)};
    coco_plot_bd(thmi,rnames1{i},'omega','w');
end
xlabel('$\tilde{\omega}$'); ylabel('$\tilde{w}$'); hold off
lnames = sprintfc('$\\sigma$=%0.1f',sigma);
legend(['mean stab.' lnames]); axis([0.1 2 0 0.8])
annotation(gcf, 'textbox', [0.0, 0.92, 0.05, 0.05], ...
    'String', '(a)','LineStyle', 'none','Interpreter','latex');
title('$\zeta = 0.05$')
fontsize(16,"points")
exportgraphics(gcf,'plots/turning_map_sigma.png','Resolution',400)

%%
zeta = [0.05 0.1 0.15 0.2];
rnames2 = sprintfc('sturn_zeta%i_psn1',1:4);

clrmp = orderedcolors("gem");
figure('Position',[400 400 600 400]); box on; hold on;
for i = length(zeta):-1:2
    p = patch([omA{i}, fliplr(omA{i-1})],[wA{i}, fliplr(wA{i-1})], [0 0 1]);  
    p.FaceColor = clrmp(i,:); p.HandleVisibility = 'off';
    p.FaceAlpha = 0.2; p.EdgeColor = 'none';
end
a = area(omA{1}, wA{1}); a.HandleVisibility = 'off';
a.FaceColor = clrmp(1,:); a.HandleVisibility = 'off';
a.EdgeColor = 'none'; a.FaceAlpha = 0.2;
for i = 1:length(zeta)
    % plot(omA{i}, wA{i},'Color',clrmp(i,:), 'LineWidth', 1.5);
    thmi.lspec = {'-', 'LineWidth', 1.5,'Color',clrmp(i,:)};
    coco_plot_bd(thmi,rnames2{i},'omega','w');
end
xlabel('$\tilde{\omega}$'); ylabel('$\tilde{w}$'); hold off
lnames = sprintfc('$\\zeta$=%0.3f',zeta);
annotation(gcf, 'textbox', [0.0, 0.92, 0.05, 0.05], ...
    'String', '(b)','LineStyle', 'none','Interpreter','latex');
legend(lnames); axis([0.1 2 0 1.2])
title('$\sigma = 1$')
fontsize(16,"points")
exportgraphics(gcf,'plots/turning_map_zeta.png','Resolution',400)