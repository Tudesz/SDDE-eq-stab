function zd_triangplotstab(restriang,color)
% function c=zd_triangplotstab(points,color)
% "plot the result of the triangle stability"
% INPUT:
%   restriang: result of the zd_trianglestab
%   color: color letter according to plot
% OUTPUT:
% VERSION: zd_triangplotstab 3.2
%       1.0: basic code
%       2.0: automatical setting of the resolution
%       3.0: plot the frequencies too
%       3.2: refreshed compatibility for Matlab2017
%% COPYRIGHT & DEVELOPMENT
% DATE 1: Vancouver, BC, Canada, UBC, 2007. 10. 23.
% DATE 2: Ideko, Spain (Basque Country), 2008. 11. 20.
% Zoltan Dombovari
% Budapest University of Technology and Echonomics
% Department of Applied Mechanics
% 2008
% ddombo@gmail.com
%% MAIN CODE

%select the points
for i=1:length(restriang.points)
    X(i)=restriang.points(i).x;
    Y(i)=restriang.points(i).y;
    Z(i)=restriang.points(i).res.aceig;
end
%borders
Dx=max(X)-min(X);
Dy=max(Y)-min(Y);
% calculate the mesh resolution
% base mesh resolution in y direction
Nyb=Dy/abs(restriang.points(2).y-restriang.points(1).y);
% the level of the converged triangles
lc=restriang.triangs(restriang.borders(1)).level;
% the resolution of the stability chart
Ny=Nyb*2^(lc/2);
%the required resolution of the mesh
Ng=ceil(4*Ny);
%grid;
[Xi,Yi]=meshgrid(min(X):Dx/Ng:max(X),min(Y):Dy/Ng:max(Y));
%interpolation of the stability
Zi=griddata(X,Y,Z,Xi,Yi,'linear');
%contourplot
contour(Xi,Yi,Zi,[1 1],'LineWidth',2,'EdgeColor',color);

return;
