clear; close all; clc;
%% Computational cost of second-moment stability in a linear, scalar SDDE
% dx(t) = (a*x(t)+b*x(t-tau) + c)*dt + 
%   (alpha*x(t)+beta*x(t-tau) + gamma)*dWt 
addpath(genpath('numerical methods'))


%% System definition

% Default parameter vector
par = [-3.0 1.0 0 -1.5 0.0 1 1]; % [a b c alpha beta gamma tau]
tau = par(7); d = 1;

% stochastic system definition
sysDDE.f = @(x,xd,t,p) p(1)*x(1,:) + p(2)*xd(1,:,1) + p(3); % drift term
sysDDE.g = @(x,xd,t,p) p(4)*x(1,:) + p(5)*xd(1,:,1) + p(6); % diffustion term

% Linear SDDE definition
sysLDDE.F = @(t,p) p(1);    % present drift term
sysLDDE.Fd = @(t,p) p(2);   % delayed drift term
sysLDDE.f = @(t,p) p(3);    % constant drift term
sysLDDE.G = @(t,p) p(4);    % present diffusion term
sysLDDE.Gd = @(t,p) p(5);   % delayed diffusion term
sysLDDE.g = @(t,p) p(6);    % constant diffusion term

% Discretization dimensions
Ms = [100 80 60 50 40 30 20 10 5];
% Ms = [5 10 20];
ne = 10; % number of evaluations at each step


%% Evaluation of the stationary covariance matrix

% Analytic solution 
[A, f0] = corr_stat_asol3(sysLDDE,tau,par);
phi_a = @(th) reshape([eye(d^2) zeros(d^2)]*expm(A*(th+tau/2))*[f0;f0],d,d);

a_time = NaN(1,length(Ms));
f_time = NaN(1,length(Ms));
x_time = NaN(1,length(Ms));
c_time = NaN(1,length(Ms));
err_f = NaN(1,length(Ms));
err_x = NaN(1,length(Ms));
for i = 1:length(Ms)
    theta = cheb_mesh(Ms(i),[-tau 0]);
    xt_mean = @(t) zeros(Ms(i),1);
    fprintf('\nM = %i\n',Ms(i));

    % Evaluations of the analytic solution
    rta = tic;
    for q = 1:ne
        C_psi = zeros(Ms(i),Ms(i));
        for j = 1:length(theta)
            for k = 1:length(theta)
                C_psi(j,k) = phi_a(-abs(theta(j)-theta(k)));
            end
        end
        fprintf('%i/%i: Analytic solution: %0.3f s\n',q,ne,toc(rta));
    end
    a_time(i) = toc(rta);
    
    % Evaluate eigenvalues of C
    rtc = tic;
    for j = 1:ne
        lamC = eig(C_psi);
        fprintf('%i/%i: Eigenvalues of C: %0.3f s\n',j,ne,toc(rtc));
    end
    c_time(i) = toc(rtc);
   
    
    % Reduced solution with Xi
    rtx = tic;
    for j = 1:ne
        JPh = pde_corr_jac(0,sysLDDE,par,tau,Ms(i));
        bPh = pde_corr_const(0,xt_mean,sysLDDE,par,tau,Ms(i));
        phi_x = -JPh\bPh;
        C_Xi = cov_interp(repmat(phi_x,1,Ms(i)).',theta); 
        fprintf('%i/%i: Reduced solution (Xi): %0.3f s\n',j,ne,toc(rtx));
    end
    x_time(i) = toc(rtx);
    err_x(i) = norm(C_Xi-C_psi);

    % Full solution with F
    if Ms(i) <= 100
        [Q, R, ~, ~] = vecmat_symm(Ms(i));
        rtf = tic;
        for j = 1:ne
            JC = pde_cov_jac(0,sysLDDE,par,tau,Ms(i));
            bC = pde_cov_const(0,xt_mean,sysLDDE,par,tau,Ms(i));
            C_Fv = -(Q*JC*R)\(Q*reshape(bC,[],1));
            C_F = reshape(R*C_Fv,Ms(i),Ms(i));
            fprintf('%i/%i: Full solution (F): %0.3f s\n',j,ne,toc(rtf));
        end
        f_time(i) = toc(rtf);
        err_f(i) = norm(C_F-C_psi);
    end

    % figure(); 
    % subplot(1,3,1); surf(theta,theta,C_psi); title('$\Psi$')
    % xlabel('$\theta$'); xlabel('$\vartheta$'); zlabel('$\bar{C}$')
    % subplot(1,3,2); surf(theta,theta,C_F); title('$\mathbf{F}$')
    % xlabel('$\theta$'); xlabel('$\vartheta$'); zlabel('$\bar{C}$')
    % if Ms(i) <= 100
    %     subplot(1,3,3); surf(theta,theta,C_Xi); title('$\Xi$')
    %     xlabel('$\theta$'); xlabel('$\vartheta$'); zlabel('$\bar{C}$')
    % end

end

% save data
save('data/scalar_cbar_time.mat','a_time','x_time','f_time',...
    'err_f','err_x','ne','Ms','c_time')

%% Comparison plots
load('data/scalar_cbar_time.mat')
figure(); 
loglog(Ms,err_f,Ms,err_x);
xlabel('$M$'); ylabel('err($\bar{C}$)');
legend('$\mathbf{F}$','$\Xi$')
figure(); 
loglog(Ms,a_time,Ms,f_time,Ms,x_time);
xlabel('$M$'); ylabel('CPU time (s)');
legend('$\Psi$','$\mathbf{F}$','$\Xi$')


%% Determining the non-uniqueness of solutions

% Find a point on the stability boundary
zerofun = @(a) scalar_char_sn(sysLDDE,tau,[a par(2:end)]);
a_l = fsolve(zerofun,par(1)); p_l = [a_l par(2:end)];

% Analytic solution
rta = tic;
Psi_a = sdde_char_mat(sysLDDE,tau,p_l);
a_time = toc(rta) + 0*Ms;

f_time = NaN(1,length(Ms));
x_time = NaN(1,length(Ms));
y_time = NaN(1,length(Ms));
err_f = NaN(1,length(Ms));
err_x = NaN(1,length(Ms));
err_y = NaN(1,length(Ms)); 
for i = 1:length(Ms)
    fprintf('\nM = %i\n',Ms(i));
    
    % Reduced solution with Xi
    rtx = tic;
    for j = 1:ne
        JPh = pde_corr_jac(0,sysLDDE,p_l,tau,Ms(i));
        lamP = eig(JPh);
        fprintf('%i/%i: Reduced spectrum (Xi): %0.3f s\n',j,ne,toc(rtx));
    end
    x_time(i) = toc(rtx);
    err_x(i) = min(abs(real(lamP)));

    % Full solution with F
    if Ms(i) <= 100
        [Q, R, ~, ~] = vecmat_symm(Ms(i));
        rtf = tic;
        for j = 1:ne
            JC = pde_cov_jac(0,sysLDDE,p_l,tau,Ms(i));
            lamC = eig(Q*JC*R);
            fprintf('%i/%i: Full spectrum (F): %0.3f s\n',j,ne,toc(rtf));
        end
        f_time(i) = toc(rtf);
        err_f(i) = min(abs(real(lamC)));
    end

    % Full solution with Y
    if Ms(i) <= 100
        % theta = cheb_mesh(Ms(i),[-tau 0]);
        sysL = pde_corr_dde(0,sysLDDE,p_l,tau,Ms(i));
        % sysL = pde_corr_dde2(0,sysLDDE,p_l,tau,Ms(i));
        rty = tic;
        for j = 1:ne
            % APh = dde_lin_igener(sysL,-theta(1:end-1),Ms(i));
            APh = dde_lin_igener2(sysL,tau,Ms(i));
            lamY = eig(APh);
            fprintf('%i/%i: Full spectrum (Y): %0.3f s\n',j,ne,toc(rty));
        end
        y_time(i) = toc(rty);
        err_y(i) = min(abs(real(lamY)));
    end

end

% save data
save('data/scalar_eig_time.mat','a_time','x_time','f_time','y_time',...
    'err_f','err_x','err_y','ne','Ms')

%% Comparison plots
load('data/scalar_eig_time.mat')
figure(); 
loglog(Ms,err_f,Ms,err_x,Ms,err_y);
xlabel('$M$'); ylabel('err($\bar{C}$)'); 
legend('$\mathbf{F}$','$\Xi$', '$\Upsilon$')
figure(); 
loglog(Ms,a_time,Ms,f_time,Ms,x_time,Ms,y_time);
xlabel('$M$'); ylabel('CPU time (s)');
legend('$\Psi$','$\mathbf{F}$','$\Xi$', '$\Upsilon$')