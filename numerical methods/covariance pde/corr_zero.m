function [d,F,JF] = corr_zero(prob,data,u,sys,mesh)
%CORR_ZERO COCO compatible zero function for tracking second moment 
%stability boundaries in stochastic delay differential equations
% Input:
%   prob: COCO problem structure
%   data: COCO data structure
%   u: vector of continuation variables [lambda; phi(theta); par]
%   sys: data sturcture of the linear SDDE:
%       -> F: matrix of drift terms @(t,p)
%       -> Fd: matrices of delayed drift terms {@(t,p)}
%       -> G: matrix of diffusion terms @(t,p)
%       -> Gd: matrices of delayed diffusion terms {@(t,p)}
%   mesh:
%       -> M: discretization dimension
%       -> n: system dimension
%       -> tau: function for evaluating the time delay
% Output:
%   d: copy of the input data structure data
%   F: evaluation of the zero funciton
%   JF: Jacobian of the zero function

% Unpack continuation data
lambda = u(1);
N = mesh.M*mesh.n^2;
phi = u(2:N+1);
par = u(N+2:end);
np = length(par);
tau = mesh.tau(par);

% Evaluate zero function
Jphi = pde_corr_jac(0,sys,par,tau,mesh.M);
F = [Jphi*phi - lambda*phi; phi.'*phi-1];
d = data;

% Evaluate Zero function Jacobian
Jphi_par = zeros(N,np);
for i = 1:np
    Jphi_dpi = pde_corr_jac_dp(0,sys,par,tau,mesh.M,i);
    Jphi_par(:,i) = Jphi_dpi*phi;
end
JF = [-phi Jphi Jphi_par; 0 2*phi.' zeros(1,np)];

end