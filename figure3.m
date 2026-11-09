clear; close all; clc;
%% Analytic second-moment stability maps of a linear, scalar SDDE
% dx(t) = (a*x(t)+b*x(t-tau) + c)*dt + 
%   (alpha*x(t)+beta*x(t-tau) + gamma)*dWt 
addpath(genpath('numerical methods'))

% RUN "scalar_moment_stab.m" BEFORE THIS SCIPT!

%% Analytic stability bounds

% Default parameter vector
par = [0 0 0 1 0 1 1]; % [a b c alpha beta gamma tau]
tau = par(7);       % time delay
a_lims = [-5 5];  % parameter limits of a
b_lims = [-5 5];  % parameter limits of b

% Linear SDDE definition
sysLDDE.F = @(t,p) p(1);    % present drift term
sysLDDE.Fd = @(t,p) p(2);   % delayed drift term
sysLDDE.f = @(t,p) p(3);    % constant drift term
sysLDDE.G = @(t,p) p(4);    % present diffusion term
sysLDDE.Gd = @(t,p) p(5);   % delayed diffusion term
sysLDDE.g = @(t,p) p(6);    % constant diffusion term

% Evaluation of implicit functions
n = 200; % contour plot resolution
a = linspace(a_lims(1),a_lims(2),n);
b = linspace(b_lims(1),b_lims(2),n);
ab_hopf = cosh(sqrt(a.^2-b.^2.')*tau/2)-...
    (a-b.')./sqrt(a.^2-b.^2.').*sinh(sqrt(a.^2-b.^2.')*tau/2);
ab_hopf(isnan(ab_hopf)) = 0;

% Do a sweep in different values of alpha
pars = [0 0; 1 0; 2 0; 1 -2; 1 -1; 1 1];
Ceq = cell(1,6);
for k = 1:6
    Ceq{k} = zeros(length(b),length(a));
    for i = 1:size(Ceq{k},1)
        for j = 1:size(Ceq{k},2)
            p_temp = par; p_temp(1:2) = [a(j) b(i)];
            p_temp(4:5) = pars(k,:);
            Ceq{k}(i,j) = scalar_char_sn(sysLDDE,tau,p_temp);
        end
    end
end


%% Plot results
load('data/scalar_stab_map.mat','rest_m')
lab = {'(a)','(b)','(c)','(d)','(e)','(f)'};

for i = 1:6
    figure('Position',[400 200 500 400]); 
    plot_stab_lambda(rest_m,1); hold on
    h = findobj(gca, 'Type', 'Contour');
    set(h,'EdgeColor','none','HandleVisibility','off')
    patch(NaN, NaN, [0.8 0.2 0.2],'FaceColor','black','EdgeColor','none',...
        'FaceAlpha',0.2);
    plot(a,a,'k:',a,-a,'k:','HandleVisibility','off')
    plot(a,-b,'k','LineWidth',1.5,'HandleVisibility','off');
    contour(a,b,ab_hopf,[0 0],'k','LineWidth',1.5,'HandleVisibility','off'); 
    contour(a,b,Ceq{i},[0 0],'LineColor','r','LineWidth',1.5); 
    [~, h1] = contourf(a,b,Ceq{i},[0 inf],'EdgeColor','none',...
        'HandleVisibility','off');
    axis([a_lims b_lims]); xlabel('$a\tau$'); ylabel('$b\tau$'); 
    annotation(gcf, 'textbox', [0.0, 0.9, 0.05, 0.05], ...
        'String', lab(i),'LineStyle', 'none','Interpreter','latex')
    title(sprintf(['$\\alpha=%0.0f \\sqrt{\\tau}$,' ...
        ' $\\beta=%0.0f \\sqrt{\\tau}$'],pars(i,1),pars(i,2))); 
    box on; fontsize(20,"points");
    drawnow
    face1 = h1.FacePrims(1);
    face1.ColorType = 'truecoloralpha';
    face1.ColorData = uint8([255; 0; 0; 20]); 
    fname = sprintf('plots/scalar_map_anal%i.png',i);
    exportgraphics(gcf,fname,'Resolution',400)
end


 

