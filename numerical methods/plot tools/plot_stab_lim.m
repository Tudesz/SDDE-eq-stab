function plot_stab_lim(restriang,m_type,ls)
%PLOT_STAB_LIM Plot stability limit based on the output of a fractal search
% -- assuming res.aceig stores charachteristic exponents --
% Input:
%   restriang: output data structure from fractal method
%   m_type: stability map type, res.aceig stores charachteristic:
%       1) exponents (limit at 0)
%       2) multipliers (limit at 1)
%   ls: linestyle data forwarded to the plot function (default empty)

if nargin<3 
    ls = {}; % no linestyle data
elseif ~iscell(ls)
    ls = {ls};
end

%select the points
X = zeros(1,length(restriang.points)); Y = X; Z = X;
for i=1:length(restriang.points)
    X(i) = restriang.points(i).x;
    Y(i) = restriang.points(i).y;
    Z(i) = restriang.points(i).res.aceig;
end

% grid with minimal required resolution
dx_min = min(diff(unique(X,'sorted')));
dy_min = min(diff(unique(Y,'sorted')));
xi = min(X):dx_min/2:max(X);
yi = min(Y):dy_min/2:max(Y);
[Xi,Yi] = meshgrid(xi,yi);

%interpolation of the stability bounds
xy = (max(X)-min(X))/(max(Y)-min(Y));
Zi = griddata(X,xy*Y,Z,Xi,xy*Yi,'linear');
if m_type == 1
    contour(Xi,Yi,Zi,[0 0],ls{:});
else
    contour(Xi,Yi,Zi,[1 1],ls{:});
end
xlim([min(X),max(X)]);

end

