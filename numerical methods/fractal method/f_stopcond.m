function res=f_stopcond(triang,points,xyrect,par)
% function res=f_stopcond(triang,points,xyrect,par)
% it is a defined function for triangle stability calculation, it gives the condition to jump out from the calculation
% INPUT:
%   triang: actual triangle,
%   points: involved points,
%   xyrect: original resolution,
%   par: system paraneters
% OUTPUT:
%   res: (0: stop, 1: separate)
% VERSION: f_stopcond 2.0
%       1.0: basic code
%       1.1: no semicolon
%       2.0: parameters included
%
% DATE 1: Vancouver, BC, Canada, UBC, 2007. 10. 23.
% DATE 2: Ideko, Spain (Basque Country), 2008. 09. 05.
% DATE 3: IK4-Ideko, Spain (Basque Country), 2017
% Zoltan Dombovari
% Budapest University of Technology and Echonomics
% Department of Applied Mechanics
% 2008
% ddombo@gmail.com
%% MAIN CODE
% global Ts;
% the necessary level is, 
% if a0 and al is the original resolution and the desired resolution

% level=2*log2(a0/al)

% tic;
    if triang.level>=5
        res=0;
    else
        res=1; 
    end
return