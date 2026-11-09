function JC = pde_cov_jac(t,sys,p0,tau,M)
%PDE_COV_JAC Evaluate the Jacobian of the second moment equation that
%describes the linear stochastic DDE
% dx(t) = (Ft*x(t) + Fdt*x(t-tau) + f)*dt + 
% (Gt*x(t) + Gdt*x(t-tau) + gt)*dWt
% Input:
%   t: integration time
%   sys: data sturcture of the linear SDDE:
%    -> F: matrix of drift terms @(t,p)
%    -> Fd: matrices of delayed drift terms {@(t,p)}
%    -> G: matrix of diffusion terms @(t,p)
%    -> Gd: matrices of delayed diffusion terms {@(t,p)}
%   p0: system parameter vector
%   tau: vector of time delays
%   M: dimension of the Chebyshev mesh used for discretization
% Output:
%   JC: Jacobian of the second moment equation wrt C

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

% Differenatiation of the moment matrix
Dm = 1/tau_max*kron(cheb_diff(M),eye(d));
A = Dm; A(i0,:) = Fm; 
B = zeros(d*M); B(i0,:) = Gm;
JC = kron(eye(M*d),A) + kron(A,eye(M*d)) +... % F*C + C*F.'
    kron(eye(M*d),B)*kron(B,eye(M*d)); % G*C*G.'

end