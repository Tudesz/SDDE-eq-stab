function dy = corr_stat_dde(t,y,Z,sys,p0,x0)
%CORR_STAT_DDE Find the derivative of the stationary correlation function
%if theta < -tau
% Input:
%   t: current time
%   y: inverted representation of the correlation function 
%       y(t) = phi(-t-tau)
%   Z: delayed evaluations of y (advanced of phi)
%   sys: data sturcture of the linear SDDE:
%    -> F: matrix of drift terms @(t,p)
%    -> Fd: matrices of delayed drift terms {@(t,p)}
%    -> f: vector of inhomogenous drift terms @(t,p)
%    -> G: matrix of diffusion terms @(t,p)
%    -> Gd: matrices of delayed diffusion terms {@(t,p)}
%    -> g: vector of inhomogenous diffusion terms @(t,p)
%   p0: system parameter vector
%   x0: stationary mean solution
% Output:
%   dy: time derivative of y at t

if ~iscell(sys.Fd)
    sys.Fd = {sys.Fd};
end
Id = eye(size(sys.F(t,p0),1));

% Unpack function definitions
dy = kron(Id,sys.F(t,p0))*y+kron(Id,sys.f(t,p0))*x0;
for i = 1:length(sys.Fd)
    dy = dy + kron(Id,sys.Fd{i}(t,p0))*Z(:,i);
end

end