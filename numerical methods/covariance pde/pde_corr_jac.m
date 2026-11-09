function JPh = pde_corr_jac(t,sys,p0,tau,M)
%PDE_CORR_JAC Evaluate the Jacobian of the correlation function equation 
% that describes the linear stochastic DDE
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
%   JPh: Jacobian of the second moment equation wrt C

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
Dm = 1/tau_max*cheb_diff(M); Dm(end,:) = 0;
JPh = kron(Dm,eye(d^2)); % partial derivatives wrt theta
theta = cheb_mesh(M,[-tau_max,0]);
Ptrans = vecmat_transp(d);

% Fill up the Jacobian with PDE rhs
for i = 1:M
    ii = (i-1)*d^2+1:i*d^2; % indicies of phi(theta)
    JPh(ii,ii) = JPh(ii,ii) + kron(eye(d),Ft); % a*phi(theta)
    for j = 1:nt
        if theta(i) + tau(j) < 0
            tqij = (theta(i)+tau(j))/tau_max + 1;
            Fij = kron(eye(d),Fdt{j}); % b*phi(theta+tau)
        else
            tqij = (-theta(i)-tau(j))/tau_max + 1;
            Fij = kron(eye(d),Fdt{j})*Ptrans; % b*phi^T(-theta-tau)
        end
        fij = lag_coeff(M,tqij);
        JPh(ii,:) = JPh(ii,:) + kron(fij,Fij);
    end
end

% Additional deterministic terms at the boundary theta = 0
i0 = (M-1)*d^2+1:M*d^2;
JPh(i0,i0) = JPh(i0,i0) + kron(Ft,eye(d)); % phi(0)*a^T
for i = 1:nt
    tqi = -tau(i)/tau_max + 1;
    fi = lag_coeff(M,tqi);
    JPh(i0,:) = JPh(i0,:) + kron(fi,kron(Fdt{i},eye(d))); % phi(-tau)*b^T
end

% Additional noise induced terms at the boundary theta = 0
tau_e = [0 tau]; Gte = [Gt, Gdt];
for i = 1:nt+1
    for j = 1:nt+1
        if tau_e(i) < tau_e(j)
            tqij = (tau_e(i)-tau_e(j))/tau_max + 1; % phi(tau_1-tau_2)
            Gij = kron(eye(d),Gte{i})*kron(Gte{j},eye(d));
        else
            tqij = (tau_e(j)-tau_e(i))/tau_max + 1; % phi^T(tau_2-tau_1)
            Gij = kron(eye(d),Gte{i})*kron(Gte{j},eye(d))*Ptrans;
        end
        fij = lag_coeff(M,tqij);
        JPh(i0,:) = JPh(i0,:) + kron(fij,Gij);
    end
end

% Avoid redundant equations
[Q, ~, iu, il] = vecmat_symm(d); nr = d^2-size(Q,1);
JPh(i0,:) = [Q*JPh(i0,:); zeros(nr,M*d^2)];
JPh(end-nr+1:end,i0(iu)) = 1;
JPh(end-nr+1:end,i0(il)) = -1;

end