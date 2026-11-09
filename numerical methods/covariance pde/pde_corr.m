function dy = pde_corr(t,y,yd,x0,sys,p0,tau)
%PDE_CORR Discretization for the froward simulation of a PDE that
%describes the correlation function of a linear DDE 
% dx(t) = (Ft*x(t) + Fdt*x(t-tau) + f)*dt + 
% (Gt*x(t) + Gdt*x(t-tau) + gt)*dWt
% Input:
%   t: integration time
%   y: current evaluation of the covariance function phi_t(theta) d*M x d
%   yd: delayed evaluation of the correlation function phi_{t-tM}
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
Ft = sys.F(t,p0);
Gt = sys.G(t,p0);
if ~iscell(sys.Fd)
    Fdt = {sys.Fd(t,p0)};
else
    Fdt = cell(1,length(tau));
    for i = 1:length(tau)
        Fdt{i} = sys.Fd{i}(t,p0);
    end
end
if ~iscell(sys.Gd)
    Gdt = {sys.Gd(t,p0)};
else
    Gdt = cell(1,length(tau));
    for i = 1:length(tau)
        Gdt{i} = sys.Gd{i}(t,p0);
    end
end
tau_max = max(tau);

% Initialization
d = size(Ft,1);
M = size(y,1)/d^2;
Pht = reshape(y,d,M*d);
Phdt = reshape(permute(reshape(yd,d,d*M,M-1),[1,3,2]),(M-1)*d,M*d);
i0 = (M-1)*d+1:M*d;

% DDE RHS evaluations
Fm = zeros(d,M*d); Fm(:,i0) = Ft;
Gm = zeros(d,M*d); Gm(:,i0) = Gt;
for i = 1:length(tau)
    tqi = -tau(i)/tau_max + 1;
    fi = lag_coeff(M,tqi);
    Fm = Fm + kron(fi,Fdt{i});
    Gm = Gm + kron(fi,Gdt{i});
end
xt = x0(t);
ft = sys.f(t,p0);
gt = sys.g(t,p0);

% Reproduction of the covariance matrix
theta = cheb_mesh(M,[-tau_max,0]);
Ct = cov_interp([Phdt; Pht],theta);

% Differenatiation of the correlation function
Dm = 1/tau_max*kron(cheb_diff(M),eye(d));
dPh = Fm*Ct + Pht*Dm.' + ft*xt.';
dPh(:,i0) = Fm*Pht.' + Pht*Fm.' + ft*xt(i0).' + xt(i0)*ft.' + ...
    Gm*Ct*Gm.' + Gm*xt*gt.' +  gt*xt.'*Gm.' + gt*gt.';
dy = reshape(dPh,[],1);

end