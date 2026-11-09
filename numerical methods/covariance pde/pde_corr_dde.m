function sysL = pde_corr_dde(t,sys,p0,tau,M)
%PDE_CORR_DDE Evaluate the Jacobian of the correlation function equation 
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
%   sysL: data sturcture of the linear DDE:
%    -> F: matrix of drift terms
%    -> Fd: matrices of delayed drift terms
%    -> f: vector of inhomogenous drift terms

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
theta = cheb_mesh(M,[-tau_max,0]);
Ptrans = vecmat_transp(d);

% Non-delayed terms (Lambda)
sysL.F = kron(Dm,eye(d^2)) +  kron(eye(d*M),Ft); % D*I + I*a
i0 = (M-1)*d^2+1:M*d^2;
sysL.F(i0,i0) = sysL.F(i0,i0) + kron(Ft,eye(d)); % a*I
sysL.F(i0,i0) = sysL.F(i0,i0) + kron(Gt,Gt); % alpha*alpha
for i = 1:nt
    tqi = -tau(i)/tau_max + 1;
    fi = lag_coeff(M,tqi);
    sysL.F(i0,:) = sysL.F(i0,:) + kron(fi,kron(Fdt{i},eye(d))) ...
        + kron(fi,kron(eye(d),Fdt{i})*Ptrans); % b*I + I*b*P
    sysL.F(i0,:) = sysL.F(i0,:) + kron(fi,kron(Gdt{i},Gt)) ...
        + kron(fi,kron(Gt,Gdt{i})*Ptrans); % beta*alpha + alpha*beta*P
end

% Additional noise induced terms at the boundary theta = 0
LFd = zeros(d^2*M,d^2*M,M);
for i = 1:M-1
    ii = (i-1)*d^2+1:i*d^2; % indicies of phi(theta)
    for j = 1:nt
        if theta(i) + tau(j) < 0
            tqij = (theta(i)+tau(j))/tau_max + 1;
            Fij = kron(eye(d),Fdt{j}); % b*phi(theta+tau)
        else
            tqij = (-theta(i)-tau(j))/tau_max + 1;
            Fij = kron(eye(d),Fdt{j})*Ptrans; % b*phi^T(-theta-tau)
        end
        fij = lag_coeff(M,tqij);
        fij2 = lag_coeff(M,1-tqij);
        LFd(ii,:,:) = LFd(ii,:,:) + kron(fij,Fij) .* reshape(fij2,1,1,M);
    end
end
for i = 1:nt
    tqi = -tau(i)/tau_max + 1;
    fi = lag_coeff(M,tqi);
    LFd(i0,i0,:) = LFd(i0,i0,:) + kron(Gdt{i},Gdt{i}) .* reshape(fi,1,1,M);
end
sysL.Fd = squeeze(num2cell(LFd(:,:,1:end-1), [1,2]));

end