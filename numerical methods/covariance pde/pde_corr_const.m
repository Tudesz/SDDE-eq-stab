function bPh = pde_corr_const(t,x0,sys,p0,tau,M)
%PDE_CORR_CONST Evaluate the constant part of the correlation function
% equation that describes the linear stochastic DDE
% dx(t) = (Ft*x(t) + Fdt*x(t-tau) + f)*dt + 
% (Gt*x(t) + Gdt*x(t-tau) + gt)*dWt
% Input:
%   t: integration time
%   x0: mean solution as a function of t @(t)
%   sys: data sturcture of the linear SDDE:
%    -> F: matrix of drift terms @(t,p)
%    -> Fd: matrices of delayed drift terms {@(t,p)}
%    -> G: matrix of diffusion terms @(t,p)
%    -> Gd: matrices of delayed diffusion terms {@(t,p)}
%   p0: system parameter vector
%   tau: vector of time delays
%   M: dimension of the Chebyshev mesh used for discretization
% Output:
%   bPh: constant part of the second moment equation

% Unpack function definitions
nt = length(tau);
Gt = sys.G(t,p0);
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
d = size(Gt,1);
i0 = (M-1)*d+1:M*d;

% DDE RHS evaluations
Gm = zeros(d,M*d); Gm(:,i0) = Gt;
for i = 1:nt
    tqi = -tau(i)/tau_max + 1;
    fi = lag_coeff(M,tqi);
    Gm = Gm + kron(fi,Gdt{i});
end
xt = x0(t);
ft = sys.f(t,p0);
gt = sys.g(t,p0);

% Differenatiation of the moment matrix
bPh = xt*ft.';
bPh(i0,:) = ft*xt(i0).' + xt(i0)*ft.' +...
    + Gm*xt*gt.' +  gt*xt.'*Gm.' + gt*gt.';
bPh = reshape(bPh,[],1);

% Avoid redundant equations
i0 = (M-1)*d^2+1:M*d^2;
Q = vecmat_symm(d);
bPh(i0,:) = [Q*bPh(i0,:); zeros(d^2-size(Q,1),1)];

end