function ph = corr_interp(C,theta,t)
%CORR_INTERP Recreate a correaltion function of a stochastic DDE via 
% interpolation on the corresponding covariance operator
% Input:
%   C: covariance matrix at t: E[x(t+theta1)*x^T(t+theta2)]
%   theta: Chebyshev mesh points on [-tau_max, 0]
%   t: time querry point, must be from [-tau_max, 0]
% Output:
%   ph: @(t) correlation functions at t E[x(t)*x^T(t+theta)]

% Initialization
M = length(theta); 
tau_max = abs(theta(1));
n = size(C,1)/M;
ph = zeros(n*M,n);

% Find C(t,theta2)
tq = 1+t/tau_max;
lt = lag_coeff(M,tq);
Ct = C*kron(lt,eye(n)).';

% Fill up the correlation function
i_min = find(t+theta<=-tau_max,1,'last');
for i = i_min:M
    ii = (i-1)*n+(1:n);
    tq = 1+(theta(i)+t)/tau_max;
    lij = lag_coeff(M,tq);
    ph(ii,:) = kron(lij,eye(n))*Ct;
end


end