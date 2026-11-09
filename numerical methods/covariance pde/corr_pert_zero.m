function [d,F] = corr_pert_zero(prob,data,u,sys,tau,M,n)
%CORR_PERT_ZERO COCO compatible zero function for tracking second moment 
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
%   tau: vector of time delays
%   M: discretization dimension
%   n: system dimension
% Output:
%   d: copy of the input data structure data
%   F: evaluation of the zero funciton

% Unpack continuation data
lambda = u(1);
phi = u(2:1+M*n^2);
par = u(M*n^2+2:end);

% Evaluate zero function
Jc = corr_pert_jac(lambda,sys,tau,par,M);
F = Jc*phi;
d = data;


end