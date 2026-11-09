function dy = pde_cov(t,y,x0,sys,p0,tau)
%PDE_COV Discretization for the froward simulation of a PDE that
%describes the second moment "function" of a linear DDE 
% dx(t) = (Ft*x(t) + Fdt*x(t-tau) + f)*dt + 
% (Gt*x(t) + Gdt*x(t-tau) + gt)*dWt
% Input:
%   t: integration time
%   y: current evaluation of the covariance function C_t(t1,t2) d*M x d*M
%   x0: mean solution as a function of t @(t)
%   sys: data sturcture of the linear SDDE:
%    -> F: matrix of drift terms @(t,p)
%    -> Fd: matrices of delayed drift terms {@(t,p)}
%    -> f: vector of inhomogenous drift terms @(t,p)
%    -> G: matrix of diffusion terms @(t,p)
%    -> Gd: matrices of delayed diffusion terms {@(t,p)}
%    -> g: vector of inhomogenous diffusion terms @(t,p)
%   p0: system parameter vector
%   tau: vector of time delays
% Output:
%   dy: time derivative of C_t(t1,t2)

% Unpack function definitions
nt = length(tau);
Ft = sys.F(t,p0);
Gt = sys.G(t,p0);
if ~iscell(sys.Fd)
    Fdt = {sys.Fd(t,p0)};
else
    Fdt = cell(1,nt);
    for i = 1:nt
        Fdt{i} = sys.Fd{i}(t,p0);
    end
end
if ~iscell(sys.Gd)
    Gdt = {sys.Gd(t,p0)};
else
    Gdt = cell(1,nt);
    for i = 1:nt
        Gdt{i} = sys.Gd{i}(t,p0);
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
Gm = zeros(d,M*d); Gm(:,i0) = Gt;
for i = 1:nt
    tqi = -tau(i)/tau_max + 1;
    fi = lag_coeff(M,tqi);
    Fm = Fm + kron(fi,Fdt{i});
    Gm = Gm + kron(fi,Gdt{i});
end
xt = x0(t);
ft = sys.f(t,p0);
gt = sys.g(t,p0);

% Differenatiation of the moment matrix
Dm = 1/tau_max*kron(cheb_diff(M),eye(d));
dCt = Dm*Ct + Ct*Dm.';
dCt(i0,:) = Fm*Ct + Ct(i0,:)*Dm.' + ft*xt.';
dCt(:,i0) = Ct*Fm.'+ Dm*Ct(:,i0) + xt*ft.';
dCt(i0,i0) = Fm*Ct(:,i0) + Ct(i0,:)*Fm.' + ft*xt(i0).' + xt(i0)*ft.' +...
    Gm*Ct*Gm.' + Gm*xt*gt.' +  gt*xt.'*Gm.' + gt*gt.';
dy = reshape(dCt,[],1);

end