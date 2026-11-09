function [d,F,JF] = cov_red_zero(prob,data,u,sys,mesh)
%COV_RED_ZERO COCO compatible zero function for tracking second moment 
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
%       -> Q: right symmetry projection matrix
%       -> R: left symmetry projection matrix
% Output:
%   d: copy of the input data structure data
%   F: evaluation of the zero funciton
%   JF: Jacobian of the zero function

% Unpack continuation data
lambda = u(1);
% N = mesh.M^2*mesh.n^2;
N = size(mesh.Q,1);
phi = u(2:N+1);
par = u(N+2:end);
np = length(par);
tau = mesh.tau(par);

% Evaluate zero function
% JC = pde_cov_jac(0,sys,par,tau,mesh.M);
JC = mesh.Q*pde_cov_jac(0,sys,par,tau,mesh.M)*mesh.R;
F = [JC*phi - lambda*phi; phi.'*phi-1];
d = data;

% Evaluate Zero function Jacobian
JC_par = zeros(N,np);
for i = 1:np
    % JC_dpi = pde_cov_jac_dp(0,sys,par,tau,mesh.M,i);
    JC_dpi = mesh.Q*pde_cov_jac_dp(0,sys,par,tau,mesh.M,i)*mesh.R;
    JC_par(:,i) = JC_dpi*phi;
end
JF = [-phi JC JC_par; 0 2*phi.' zeros(1,np)];

end