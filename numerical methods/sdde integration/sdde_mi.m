function [t,U,Ud] = sdde_mi(f,g,dg,p,U0,tau,tw,dW,opts)
%SDDE_MI Milstein integration routine for stochastic DDEs (evaluating
%an ensamble of solutions in paralell)
%TODO: Extend to non-scalar systems
% Input:
%   f: ode rhs diffusion term @(x,xd,t,p) 
%   g: ode rhs drift terms, cell array: {@(x,xd,t,p)}
%   dg: ode rhs drift Jacobians wrt x, cell array: {@(x,xd,t,p)}
%   p: system parameter vector
%   U0: initial distribution sample of [x(t0) std(x(t0))]
%   tau: fixed, point time delays
%   tw: time signal of the brownian paths
%   dW: increments of the a priori simulated Brownian paths
%   opts: options data structure
%    -> R: number of Brownian increments between Euler steps (default 1)
%    -> logs: if true log the progress of the integration (default false)
% Output:
%   t: integration time mesh (length(tw)/opts.R)
%   U: solution realisations evaluated on t

if nargin < 8 || ~isfield(opts,'R')
    R = 1;
else
    R = opts.R;
end
if nargin < 8 || ~isfield(opts,'logs')
    opts.logs = false;
end

% Output initialization
t = tw(1:R:end);
N = length(t);
d = size(U0(0),1);
n = size(dW,1);
U = zeros(d,N,n);
m = length(g);
nt = length(tau);
Ud = zeros(d,N,n,nt);

% Forward time integration
if opts.logs
    fprintf('\nSDDE solution using a fix step Milstein method\n')
    tic;
end
U(:,1,:) = U0(repmat(t(1),1,n));
dt = tw(end)/(length(tw)-1);
DT = R*dt; % time increment
for i = 1:N

    % find the delayed terms
    xdi = zeros(d,n,nt);
    for j = 1:nt
        tdel = t(i) - tau(j);
        if tdel > t(1)
            id = i-ceil(tau(j)/DT);
            xdi(:,:,j) = reshape(U(:,id,:),d,[]);
        else
            xdi(:,:,j) = reshape(U0(repmat(tdel,1,n)),d,[]);
        end
    end
    Ud(:,i,:,:) = xdi;

    if i>1
        % do an Euler-Maruyama step
        dW_ind = 1+R*(i-2) + (1:R);
        xi = reshape(U(:,i-1,:),d,[]);
        F = f(xi,xdi,t(i),p); % drift term
        for j = 1:m
            G = g{j}(xi,xdi,t(i),p); % diffusion term
            JG = dg{j}(xi,xdi,t(i),p); % diffusion term Jacobian
            DW = reshape(sum(dW(:,dW_ind,j),2),1,[]); % Wiener increment
            xi = xi + G.*DW + 1/2*G.*JG.*(DW.^2-DT); 
        end
        U(:,i,:) = xi + F.*DT;
        if opts.logs && mod(i,floor(N/10)) == 0 && i > 1
            fprintf('   solution progress: %0.1f%%, |X| = %d, time: %0.3f s\n',...
                i/N*100,norm(squeeze(U(:,i,:))),toc);
        end
    end
end