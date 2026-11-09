function [t, W, dW] = brownian_path(T,dt,n,m)
%BROWNIAN_PATH Simulate a trajectory of the standard Wiener process
% Input:
%   T: signal length
%   dt: timestep
%   n: number of trajectories generated
%   m: number of non-correlated path needed (default 1)
% Output:
%   t: time signal
%   W: Wiener process on t
%   dW: increments of the Wiener process

if nargin < 4
    m = 1;
end

% Initialize signals
t = 0:dt:T;
N = length(t);

% Evaluate the random variables
dW = cat(2,zeros(n,1,m),sqrt(dt)*randn(n,N-1,m));
W = cumsum(dW,2);

end

