function sysL = pde_corr_dde2(t,sys,p0,tau,M)
%PDE_CORR_DDE2 Evaluate the Jacobian of the correlation function equation 
% that describes the linear stochastic DDE (only works with one delay!)
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
Ft = sys.F(t,p0);
Gt = sys.G(t,p0);
if ~iscell(sys.Fd)
    Fdt = sys.Fd(t,p0);
else
    Fdt = sys.Fd{1}(t,p0);
end
if ~iscell(sys.Gd)
    Gdt = sys.Gd(t,p0);
else
    Gdt = sys.Gd{1}(t,p0);
end

% Initialization
d = size(Ft,1);
Dm = 1/tau*cheb_diff(M); Dm(end,:) = 0;
Ptrans = vecmat_transp(d);

% Non-delayed terms (Lambda)
sysL.F = kron(Dm,eye(d^2)) +  kron(eye(d*M),Ft); % D*I + I*a
i0 = (M-1)*d^2+1:M*d^2;
sysL.F(i0,i0) = sysL.F(i0,i0) + kron(Ft,eye(d)); % a*I
sysL.F(i0,i0) = sysL.F(i0,i0) + kron(Gt,Gt); % alpha*alpha
fi = zeros(1,M); fi(1) = 1;
sysL.F(i0,:) = sysL.F(i0,:) + kron(fi,kron(Fdt,eye(d))) ...
    + kron(fi,kron(eye(d),Fdt)*Ptrans); % b*I + I*b*P
sysL.F(i0,:) = sysL.F(i0,:) + kron(fi,kron(Gdt,Gt)) ...
    + kron(fi,kron(Gt,Gdt)*Ptrans); % beta*alpha + alpha*beta*P


% Additional noise induced terms at the boundary theta = 0
sysL.Fd = cell(1,M-1);
for i = 1:M-1
    sysL.Fd{i} = zeros(d^2*M,d^2*M);
    ii = (i-1)*d^2+1:i*d^2; % indicies of phi(theta)
    Fi = kron(eye(d),Fdt)*Ptrans; % b*phi^T(-theta-tau)
    fi = zeros(1,M); fi(M+1-i) = 1;
    sysL.Fd{i}(ii,:) = sysL.Fd{i}(ii,:) + kron(fi,Fi);
end
sysL.Fd{1}(i0,i0) = sysL.Fd{1}(i0,i0) + kron(Gdt,Gdt);

end