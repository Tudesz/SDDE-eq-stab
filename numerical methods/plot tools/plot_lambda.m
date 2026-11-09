function plot_lambda(lam)
%PLOT_LAMBDA Plot characteristic eigenvalues on the Im-Re plane
% Input:
%   lam: array of characteristic multipliters

scatter(real(lam),imag(lam),50,'bx'); hold on
scatter(real(lam(real(lam)>0)),imag(lam(real(lam)>0)),50,'rx'); 
xline(0,':k'); hold off
xlabel('$\Re$'); ylabel('$\Im$'); 
[Re_max,ind] = max(real(lam));
title(sprintf('Critical eigenvalue: %0.5f $\\pm$ %0.5fi',...
    real(lam(ind)),imag(lam(ind))))
xlim([-inf  max([0.1*max(abs(real(lam))),Re_max])]);
box on;
[~,ic] = max(real(lam));
title(sprintf('$\\lambda_{\\mathrm{crit}}$ = %0.3f + %0.3f i',...
    real(lam(ic)),imag(lam(ic))));
end

