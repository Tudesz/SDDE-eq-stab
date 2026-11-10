clear; close all; clc;
%% Computational cost of second-moment stability in a linear, scalar SDDE
% dx(t) = (a*x(t)+b*x(t-tau) + c)*dt + 
%   (alpha*x(t)+beta*x(t-tau) + gamma)*dWt 
addpath(genpath('numerical methods'))

% RUN "calculations/scalar_comp_time.m" BEFORE THIS SCIPT!


%% Preprocess convergence study data

% Evaluation of C
load('data/scalar_cbar_time.mat','a_time','x_time','f_time',...
    'err_f','err_x','ne','Ms','c_time')
Cat = a_time./ne; Cft = f_time./ne; Cxt = x_time./ne; 
Cfe = err_f; Cxe = err_x;
Catp = polyfit(log(Ms),log(Cat),1);
Cftp = polyfit(log(Ms),log(Cft),1); 
Cfep = polyfit(log(Ms),log(Cfe),1);
Cxtp = polyfit(log(Ms),log(Cxt),1); 
Cxep = polyfit(log(Ms(end-2:end)),log(Cxe(end-2:end)),1);
Ct = c_time./ne; Ctp = polyfit(log(Ms),log(Ct),1);
fprintf('\neig(C) ~ M^{%0.3f}\n\n',Ctp(1))

load('data/scalar_eig_time.mat','a_time','x_time','f_time','y_time',...
    'err_f','err_x','err_y','ne','Ms')
Fat = a_time./ne; Fft = f_time./ne; Fxt = x_time./ne; Fyt = y_time./ne;
Ffe = err_f; Fxe = err_x; Fye = err_y;
Fatp = polyfit(log(Ms),log(Fat),1);
Fftp = polyfit(log(Ms),log(Fft),1); 
Ffep = polyfit(log(Ms),log(Ffe),1);
Fxtp = polyfit(log(Ms),log(Fxt),1); 
Fxep = polyfit(log(Ms(end-1:end)),log(Fxe(end-1:end)),1);
Fyep = polyfit(log(Ms(end-1:end)),log(Fye(end-1:end)),1);
Fytp = polyfit(log(Ms),log(Fyt),1); 

T = table([0; Cfep(1); Cxep(1)],[Catp(1); Cftp(1); Cxtp(1)],...
    [0; Ffep(1); Fxep(1)],[Fatp(1); Fftp(1); Fxtp(1)],...
    'VariableNames',{'C_err','C_time','F_err','F_time'},...
    'RowNames',{'Psi','F','Xi'});
disp(T)


%% Combined plots
figure('Position',[400 200 600 750]); 

ax1 = subplot(3,2,1);
loglog(Ms,Ffe,'bo-',Ms,Fye,'ro-','LineWidth',1.5); hold on
% loglog(Ms,exp(polyval(Ffep,log(Ms))),':b','LineWidth',1.5);
% loglog(Ms(end-1:end),exp(polyval(Fyep,log(Ms(end-1:end)))),':r','LineWidth',1.5);
xlabel('$M$'); ylabel('$|\rho(\hat{\mathbf{F}},\mathbf{\Upsilon})|$'); grid on; hold off
title([sprintf('$\\hat{\\mathbf{F}}\\sim M^{%0.1f}$',Ffep(1)), ', '...
    sprintf('$\\mathbf{\\Upsilon}\\sim M^{%0.1f} $',Fyep(1))])
ax1.XTick = [10 100]; ax1.XTickLabelRotation = 0; ax1.YTick = [1e-10 1e-5];
annotation(gcf, 'textbox', [0.02, 0.92, 0.05, 0.05], ...
    'String', '(a)','LineStyle', 'none','Interpreter','latex');
ax2 = subplot(3,2,2);
loglog(Ms,Fft,'bo-',Ms,Fyt,'ro-','LineWidth',1.5); hold on
% loglog(Ms,Fat,'ko-','LineWidth',1.5);
% loglog(Ms,exp(polyval(Fftp,log(Ms))),':b','LineWidth',1.5);
% loglog(Ms,exp(polyval(Fytp,log(Ms))),':r','LineWidth',1.5);
% loglog(Ms,exp(polyval(Fatp,log(Ms))),':k','LineWidth',1.5); 
xlabel('$M$'); ylabel('CPU time (s)'); grid on; hold off
title([sprintf('$\\hat{\\mathbf{F}}\\sim M^{%0.1f}$',Fftp(1)), ', '...
    sprintf('$\\mathbf{\\Upsilon}\\sim M^{%0.1f} $',Fytp(1))])
ax2.XTick = [10 100]; ax2.XTickLabelRotation = 0;
annotation(gcf, 'textbox', [0.5, 0.92, 0.05, 0.05], ...
    'String', '(b)','LineStyle', 'none','Interpreter','latex');


ax3 = subplot(3,2,3); 
loglog(Ms,Cfe,'bo-',Ms,Cxe,'ro-','LineWidth',1.5); hold on
% loglog(Ms,exp(polyval(Cfep,log(Ms))),':b','LineWidth',1.5);
% loglog(Ms(end-2:end),exp(polyval(Cxep,log(Ms(end-2:end)))),':r','LineWidth',1.5);
xlabel('$M$'); ylabel('$||\bar{C}_M-\bar{C}||$'); grid on;  hold off
title([sprintf('$\\hat{\\mathbf{F}}\\sim M^{%0.1f}$',Cfep(1)), ', '...
    sprintf('$\\mathbf{\\Xi}\\sim M^{%0.1f} $',Cxep(1))])
ax3.XTick = [10 100]; ax3.XTickLabelRotation = 0; ax3.YTick = [1e-10 0.1];
annotation(gcf, 'textbox', [0.02, 0.62, 0.05, 0.05], ...
    'String', '(c)','LineStyle', 'none','Interpreter','latex');
ax4 = subplot(3,2,4);
loglog(Ms,Cft,'bo-',Ms,Cxt,'ro-','LineWidth',1.5); hold on
loglog(Ms,Cat,'ko-','LineWidth',1.5); 
% loglog(Ms,exp(polyval(Cftp,log(Ms))),':b','LineWidth',1.5);
% loglog(Ms,exp(polyval(Cxtp,log(Ms))),':r','LineWidth',1.5);
% loglog(Ms,exp(polyval(Catp,log(Ms))),':k','LineWidth',1.5);
xlabel('$M$'); ylabel('CPU time (s)'); grid on;  hold off
title([sprintf('$\\hat{\\mathbf{F}}\\sim M^{%0.1f}$',Cftp(1)), ', '...
    sprintf('$\\mathbf{\\Xi}\\sim M^{%0.1f} $',Cxtp(1))])
ax4.XTick = [10 100]; ax4.XTickLabelRotation = 0;
annotation(gcf, 'textbox', [0.5, 0.62, 0.05, 0.05], ...
    'String', '(d)','LineStyle', 'none','Interpreter','latex');
% legend('$\Psi$','$\hat{\mathbf{F}}$','$\Xi$')

ax5 = subplot(3,2,5);
loglog(Ms,Ffe,'bo-',Ms,Fxe,'ro-','LineWidth',1.5); hold on
% loglog(Ms,exp(polyval(Ffep,log(Ms))),':b','LineWidth',1.5);
% loglog(Ms(end-1:end),exp(polyval(Fxep,log(Ms(end-1:end)))),':r','LineWidth',1.5);
xlabel('$M$'); ylabel('$|\mathrm{eig}(\hat{\mathbf{F}},\mathbf{\Xi})|$'); grid on; hold off
title([sprintf('$\\hat{\\mathbf{F}}\\sim M^{%0.1f}$',Ffep(1)), ', '...
    sprintf('$\\mathbf{\\Xi}\\sim M^{%0.1f} $',Fxep(1))])
ax5.XTick = [10 100]; ax5.XTickLabelRotation = 0;  
annotation(gcf, 'textbox', [0.02, 0.32, 0.05, 0.05], ...
    'String', '(e)','LineStyle', 'none','Interpreter','latex');
ax6 = subplot(3,2,6);
loglog(Ms,Fft,'bo-',Ms,Fxt,'ro-','LineWidth',1.5); hold on
% loglog(Ms,Fat,'ko-','LineWidth',1.5);
% loglog(Ms,exp(polyval(Fftp,log(Ms))),':b','LineWidth',1.5);
% loglog(Ms,exp(polyval(Fxtp,log(Ms))),':r','LineWidth',1.5);
% loglog(Ms,exp(polyval(Fatp,log(Ms))),':k','LineWidth',1.5); 
xlabel('$M$'); ylabel('CPU time (s)'); grid on; hold off
title([sprintf('$\\hat{\\mathbf{F}}\\sim M^{%0.1f}$',Fftp(1)), ', '...
    sprintf('$\\mathbf{\\Xi}\\sim M^{%0.1f} $',Fxtp(1))])
ax6.XTick = [10 100]; ax6.XTickLabelRotation = 0;
annotation(gcf, 'textbox', [0.5, 0.32, 0.05, 0.05], ...
    'String', '(f)','LineStyle', 'none','Interpreter','latex');
fontsize(16,'points')
% legend('$\Psi$','$\hat{\mathbf{F}}$','$\Xi$')
% exportgraphics(gcf,'plots/scalar_convergence.png','Resolution',400)