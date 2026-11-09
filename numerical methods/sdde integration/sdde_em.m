function [t,U,Ud,dWt] = sdde_em(f,g,p,U0,tau,tw,dW,opts)
%SDDE_EM Euler-maruyama integration routine for stochastic DDEs (evaluating
%an ensamble of solutions in paralell)
% Input:
%   f: ode rhs diffusion term @(x,xd,t,p) 
%   g: ode rhs drift terms, cell array: {@(x,xd,t,p)}
%   p: system parameter vector
%   U0: initial distribution sample of [x(t0) std(x(t0))]
%   tau: fixed, point time delays
%   tw: time signal of the brownian paths
%   dW: increments of the a priori simulated Brownian paths
%   opts: options data structure
%    -> R: number of steps between saved values (default 1)
%    -> logs: if true log the progress of the integration (default false)
% Output:
%   t: integration time mesh (length(tw)/opts.R)
%   U: solution realisations evaluated on t
%   Ud: realisations of the delayed state U(t-tau)
%   dWt: realisations of the Wiener increments on t

if nargin < 7 || ~isfield(opts,'R')
    R = 1;
else
    R = opts.R;
end
if nargin < 7 || ~isfield(opts,'logs')
    opts.logs = false;
end

% Output initialization
t = tw(1:R:end);
dt = (tw(end)-tw(1))/(length(tw)-1);
N = length(t);
d = size(U0(0),1);
n = size(dW,1);
U = zeros(d,N,n);
m = length(g);
nt = length(tau);
Ud = zeros(d,N,n,nt);
dWt = zeros(N,n,m);


% Forward time integration
if opts.logs
    fprintf('\nSDDE solution using a fix step EM method\n')
    tic;
end
U(:,1,:) = U0(repmat(t(1),1,n));
xdi0 = zeros(d,n,nt);
for j = 1:nt
    xdi0(:,:,j) = reshape(U0(repmat(t(1)-tau(j),1,n)),d,[]);
end
Ud(:,1,:,:) = squeeze(xdi0);

for i = 2:N
    ind_R = R*(i-2) + (1:R);
    dW_i = dW(:,ind_R,:);
    tw_i = tw(ind_R);
    xk = reshape(U(:,i-1,:),d,[]);

    % find the delayed terms
    xdi1 = zeros(d,n,nt);
    for j = 1:nt
        tdel = t(i) - tau(j);
        if tdel > t(1)
            id = min([ceil(tdel/(R*dt)),i-1]);
            xdi1(:,:,j) = reshape(U(:,id,:),d,[]);
        else
            xdi1(:,:,j) = reshape(U0(repmat(tdel,1,n)),d,[]);
        end
    end

    for k = 1:R
        % Interpolate the delayed terms (linear)
        xdk = (R-k+1)/R*xdi0 + (k-1)/R*xdi1;

        % do an Euler-Maruyama step
        F = f(xk,xdk,tw_i(k),p); % drift term
        xk = xk + F.*dt;
        for j = 1:m
            G = g{j}(xk,xdk,tw_i(k),p); % diffusion term
            dw = squeeze(dW_i(:,k,j)).'; % Wiener increment
            xk = xk + G.*dw; 
        end

    end
    Ud(:,i,:,:) = xdi1;
    U(:,i,:) = xk;
    dWt(i,:,:) = squeeze(sum(dW_i,2));
    if opts.logs && mod(i,floor(N/10)) == 0 && i > 1
        fprintf('   solution progress: %0.1f%%, |X| = %d, time: %0.3f s\n',...
            i/N*100,norm(squeeze(U(:,i,:))),toc);
    end
    xdi0 = xdi1;
end

end