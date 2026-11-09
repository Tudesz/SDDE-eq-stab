function dy = pde_cov_proj(t,y,sys,p0,tau)
%PDE_COV Discretization for the froward simulation of a PDE that
%describes the moment "function" of a linear DDE 
% x'(t) =Ft*x(t) + Fdt*x(t-tau) + Gt*dWt
% Input:
%   t: integration time
%   y: current evaluation of the covariance function C_t(t1,t2) d*M x d*M
%   sys: data sturcture of the projected linear SDDE:
%    -> F: matrix of drift terms @(t,p)
%    -> Fd: matrices of delayed drift terms {@(t,p)}
%    -> gt: projected inhomogenous diffusion terms on theta @(t,theta,p)
%    -> g0: projected inhomogenous diffusion terms at theta = 0 @(t,p)
%   p0: system parameter vector
%   tau: vector of time delays
% Output:
%   dy: time derivative of C_t(t1,t2)

% Unpack function definitions
nt = length(tau);
Ft = sys.F(t,p0);
if ~iscell(sys.Fd)
    Fdt = {sys.Fd(t,p0)};
else
    Fdt = cell(1,nt);
    for i = 1:nt
        Fdt{i} = sys.Fd{i}(t,p0);
    end
end
tau_max = max(tau);

% Initialization
d = size(Ft,1);
M = sqrt(size(y,1)/d^2);
Ct = reshape(y,M*d,M*d);
i0 = (M-1)*d+1:M*d;

% DDE RHS evaluations
Fm = zeros(d,M*d); Fm(:,i0) = Ft;
for i = 1:nt
    tqi = -tau(i)/tau_max + 1;
    fi = lag_coeff(M,tqi);
    Fm = Fm + kron(fi,Fdt{i});
end
g0 = sys.g0(t,p0);

% Intermittent terms in the stochastic excitation
theta = cheb_mesh(M,[-tau_max 0]);
gt = zeros(d,M);
for i = 1:M
    gt(:,i) = sys.gt(t,theta(i),p0);
end
gt = reshape(gt,[],1);

% Differenatiation of the moment matrix
Dm = 1/tau_max*kron(cheb_diff(M),eye(d));
dCt = Dm*Ct + Ct*Dm.' + gt*gt.';
dCt(i0,:) = Fm*Ct + Ct(i0,:)*Dm.' + g0*gt.';
dCt(:,i0) = Ct*Fm.'+ Dm*Ct(:,i0) + gt*g0.';
dCt(i0,i0) = Fm*Ct(:,i0) + Ct(i0,:)*Fm.' + g0*g0.';
dy = reshape(dCt,[],1);

end