function A = dde_lin_igener2(sys,tau,M)
%DDE_LIN_IGENER2 Infinitesimal generator of a linear DDE (special case
% without interpolation)
% Input:
%   sys: data sturcture of the linear SDDE:
%    -> F: matrix of drift terms
%    -> Fd: matrices of delayed drift terms
%   tau: vector of time delays
%   M: Chebyshev mesh dimension
% Output:
%   A: infinitesimal generator of the local linear system


% initialize the A matrix
n = size(sys.F,1);
Dm = (M-1)/tau*(diag(ones(M-1,1),1)-eye(M));
A = kron(Dm,eye(n));

% Fill up the A matrix
A(end-n+1:end,:) = 0; % clear rows corresponding to the present
A(end-n+1:end,end-n+1:end) = sys.F; % effect of non delayed terms
for i = 1:M-1
    lq = zeros(1,M); lq(i) = 1;
    A(end-n+1:end,:) = A(end-n+1:end,:) + kron(lq,sys.Fd{i});
end


end