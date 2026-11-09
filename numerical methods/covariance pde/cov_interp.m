function C = cov_interp(ph,theta)
%COV_INTERP Recreate the covariance matrix of a stochastic DDE by 
% interpolating delayed instances of the correlation functions exploiting
% the fact that ph_t(th) = C_t(0,th) and C_t(th1,th2) = ph_{t+th1}(th2-th1)
% Input:
%   ph: matrix of correlation functions at t+tau_i: i = 1,...M
%       E[x(t+tau_i)*x^T(t+tau_i+theta)], tau_i = theta_i 
%       stored as row vectors [ph(t-tau_max); ...; ph(t)]
%   theta: Chebyshev mesh points on [-tau_max, 0]
% Output:
%   C: covariance matrix at t: E[x(t+theta1)*x^T(t+theta2)]

% Initialization
M = length(theta); 
tau_max = abs(theta(1));
n = size(ph,1)/M;

% fill up the covariance matrix considering
C = zeros(n*M); %  C_t(th1,th2)
for i = 1:M
    ii = (i-1)*n+(1:n); % index of x(th1)
    for j = 1:M
        jj = (j-1)*n+(1:n); % index of x^T(th2)
        tq = theta(i); % querry point in actual time (t+th1)
        thq = theta(j)-theta(i); % querry point in local time (th2-th1)
        if thq <= 0 % interpolation only valid if thq in [-tau_max 0]
            lt = lag_coeff(M,1+tq/tau_max); % interpolation in t
            lth = lag_coeff(M,1+thq/tau_max); % interpolation in tau
            C(ii,jj) = kron(lt,eye(n)) * ph * kron(lth,eye(n)).';
            C(jj,ii) = C(ii,jj).'; % exploit the symmetry of C_t
        end
    end
end

end