function zd_triangleplot(points,triangs,borders)
% plot the calculated triangle mesh based on the zd_trianglestab file
% INPUT:
%   points, triangs, borders: results from zd_trianglestab file
% OUTPUT: 
%% VERSION 2.0
%% COPYRIGHT & DEVELOPMENT
% Zoltan Dombovari
% Budapest University of Technology and Echonomics
% Department of Applied Mechanics
% 2007
% ddombo@gmail.com
% 2017, Ideko-IK4
%% MAIN CODE
%declaring coordinates
kx=zeros(3,1);
ky=zeros(3,1);

hold on;
% triplot(TRI,x,y);
if ~exist('borders','var')
    for j=1:length(triangs)
        switch triangs(j).type
            case 0
                sty='g';
            case 1
                sty='r';
            case 2
                sty='b';
        end
        for k=1:3
            kx(k,1)=points(triangs(j).points(k)).x;
            ky(k,1)=points(triangs(j).points(k)).y;
        end
        fill(kx,ky,sty);
%         for k=1:3
%             text(kx(k),ky(k),['-------' int2str(k)]);
%         end
%         pause;
    end
else
    sty='k';
    for j=1:length(borders)
        for k=1:3
            kx(k,1)=points(triangs(borders(j)).points(k)).x;
            ky(k,1)=points(triangs(borders(j)).points(k)).y;
        end
        fill(kx,ky,sty);
    end
end
% plot(x0,y0,'go');
% plot(x1,y1,'ro');