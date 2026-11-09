function JPh_dp = pde_corr_jac_dp(t,sys,p0,tau,M,pind)
%PDE_CORR_JAC Evaluate the parameter Jacobian of the Jacobian of the 
% correlation function equation that describes the linear stochastic DDE
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
%   JPh_dp: Jacobian of the second moment equation wrt C

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
    Gdt_p = cell(1,nt);
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
Dm = -tau_max_p/tau_max^2*cheb_diff(M); Dm(end,:) = 0;
JPh_dp = kron(Dm,eye(d^2)); % partial derivatives wrt theta
theta = cheb_mesh(M,[-tau_max,0]);
Ptrans = vecmat_transp(d);

% Fill up the Jacobian with PDE rhs
for i = 1:M
    % TODO: handle multiple delays
    ii = (i-1)*d^2+1:i*d^2; % indicies of phi(theta)
    JPh_dp(ii,ii) = JPh_dp(ii,ii) + kron(eye(d),Ft_p); % a*phi(theta)
    for j = 1:nt
        if theta(i) + tau(j) < 0
            tqij = (theta(i)+tau(j))/tau_max + 1;
            Fij = kron(eye(d),Fdt_p{j}); % b*phi(theta+tau)
        else
            tqij = (-theta(i)-tau(j))/tau_max + 1;
            Fij = kron(eye(d),Fdt_p{j})*Ptrans; % b*phi^T(-theta-tau)
        end
        fij = lag_coeff(M,tqij);
        JPh_dp(ii,:) = JPh_dp(ii,:) + kron(fij,Fij);
    end
end

% Additional deterministic terms at the boundary theta = 0
i0 = (M-1)*d^2+1:M*d^2;
JPh_dp(i0,i0) = JPh_dp(i0,i0) + kron(Ft_p,eye(d)); % phi(0)*a^T
for i = 1:nt
    tqi = -tau(i)/tau_max + 1;
    fi = lag_coeff(M,tqi);
    JPh_dp(i0,:) = JPh_dp(i0,:) + kron(fi,kron(Fdt_p{i},eye(d))); % phi(-tau)*b^T
end

% Additional noise induced terms at the boundary theta = 0
tau_e = [0 tau]; Gte_p = [Gt_p, Gdt_p]; Gte = [Gt, Gdt];
for i = 1:nt+1
    for j = 1:nt+1
        if tau_e(i) < tau_e(j)
            tqij = (tau_e(i)-tau_e(j))/tau_max + 1; % phi(tau_1-tau_2)
            Gij = kron(eye(d),Gte{i})*kron(Gte_p{j},eye(d)) + ...
                kron(eye(d),Gte_p{i})*kron(Gte{j},eye(d));
        else
            tqij = (tau_e(j)-tau_e(i))/tau_max + 1; % phi^T(tau_2-tau_1)
            Gij = kron(eye(d),Gte{i})*kron(Gte_p{j},eye(d))*Ptrans + ...
                kron(eye(d),Gte_p{i})*kron(Gte{j},eye(d))*Ptrans;
        end
        fij = lag_coeff(M,tqij);
        JPh_dp(i0,:) = JPh_dp(i0,:) + kron(fij,Gij);
    end
end

% Avoid redundant equations
[Q, ~, ~, ~] = vecmat_symm(d); nr = d^2-size(Q,1);
JPh_dp(i0,:) = [Q*JPh_dp(i0,:); zeros(nr,M*d^2)];

end