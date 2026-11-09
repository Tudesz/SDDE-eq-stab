function [d,F] = corr_pert_comp_zero(prob,data,u,sys,tau,M,n)
%CORR_PERT_ZERO COCO compatible zero function for tracking Hopf second 
%moment stability boundaries in stochastic delay differential equations
% Input:
%   prob: COCO problem structure
%   data: COCO data structure
%   u: vector of continuation variables [lambda; real(phi); imag(phi); par]
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
mu = u(1);
omega = u(2);
phi_r = u(3:2+M*n^2);
phi_i = u(3+M*n^2:2+2*M*n^2);
par = u(2*M*n^2+3:end);

% Evaluate zero function
Jc = corr_pert_jac(mu+omega*1i,sys,tau,par,M);
Jc_r = real(Jc); Jc_i = imag(Jc);
F = [Jc_r*phi_r - Jc_i*phi_i; Jc_r*phi_i + Jc_i*phi_r]; % [Re; Im]
% Fc = Jc*(phi_r+1i*phi_i);
% F = [real(Fc); imag(Fc)];
d = data;

end