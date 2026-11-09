function JC_dp = pde_cov_jac_dp(t,sys,p0,tau,M,pind)
%PDE_COV_JAC Evaluate the parameter Jacobian of the Jacobian of the second
%  moment equation that describes the linear stochastic DDE
% dx(t) = (Ft*x(t) + Fdt*x(t-tau) + f)*dt + 
% (Gt*x(t) + Gdt*x(t-tau) + gt)*dWt
% TODO: properly generlize to multiple delays!
% Input:
%   t: integration time
%   sys: data sturcture of the linear SDDE:
%    -> F_p: parameter Jacobian of the matrix of drift terms @(t,p,i)
%    -> Fd_p: parameter Jacobian of the matrices of delayed drift terms {@(t,p,i)}
%    -> G_p: parameter Jacobian of the matrix of diffusion terms @(t,p,i)
%    -> Gd_p: parameter Jacobian of the  matrices of delayed diffusion terms {@(t,p,i)}
%    -> tau_p: parameter Jacobian of the time delay vector @(t,p,i)
%   p0: system parameter vector
%   tau: vector of time delays
%   M: dimension of the Chebyshev mesh used for discretization
%   pind: index of the parameter with respect to which the Jacobain is evalutated
% Output:
%   JC_dp: parameter Jacobian of the Jacobian of the second moment equation

% Unpack function definitions
nt = length(tau);
Ft_p = sys.F_p(t,p0,pind);
Gt_p = sys.G_p(t,p0,pind);
Gt = sys.G(t,p0);
if ~iscell(sys.Fd_p)
    Fdt_p = {sys.Fd_p(t,p0,pind)};
else
    Fdt_p = cell(1,nt);
    for i = 1:nt
        Fdt_p{i} = sys.Fd_p{i}(t,p0,pind);
    end
end
if ~iscell(sys.Gd_p)
    Gdt_p = {sys.Gd_p(t,p0,pind)};
else
    Gdt_p = cell(1,nt,pind);
    for i = 1:nt
        Gdt_p{i} = sys.Gd_p{i}(t,p0,pind);
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
[tau_max, imax] = max(tau);
tau_p = sys.tau_p(t,p0,pind);
tau_max_p = tau_p(imax);

% Initialization
d = size(Ft_p,1);
i0 = (M-1)*d+1:M*d;

% DDE RHS evaluations
Fm_p = zeros(d,M*d); Fm_p(:,i0) = Ft_p;
Gm_p = zeros(d,M*d); Gm_p(:,i0) = Gt_p;
Gm = zeros(d,M*d); Gm(:,i0) = Gt;
for i = 1:nt
    % TODO: handle multiple delays
    tqi = -tau(i)/tau_max + 1;
    fi = lag_coeff(M,tqi);
    Fm_p = Fm_p + kron(fi,Fdt_p{i});
    Gm_p = Gm_p + kron(fi,Gdt_p{i});
    Gm = Gm + kron(fi,Gdt{i});
end

% Differenatiation of the moment matrix
Dm = 1/tau_max*kron(cheb_diff(M),eye(d));
A_p = -tau_max_p/tau_max*Dm;
A_p(i0,:) = Fm_p; 
B_p = zeros(d*M); B_p(i0,:) = Gm_p;
B = zeros(d*M); B(i0,:) = Gm;
JC_dp = kron(eye(M*d),A_p) + kron(A_p,eye(M*d)) +... % F*C + C*F.'
    kron(eye(M*d),B_p)*kron(B,eye(M*d)) + ...
    kron(eye(M*d),B)*kron(B_p,eye(M*d)); % G*C*G.'


end