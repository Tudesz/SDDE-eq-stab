function restriag=zd_trianglestab(condfunc,xyrect,Nbasemesh,par,stopcond,ctransf_x,ctransf_y,oldres)

% function [triangs,points,borders]=zd_trianglestab(condfunc,xyrect,Nbasemesh,stopcond)
% "this function is capable to find the stability boundaries in the given area by using
%  triangle fractals method"
% INPUT:
%   [res,stab]=condfunc(X,Y,<par>): is a function which gives back the stability (0: stable, 1: unstable)
%                    (X,Y):inner coordinates
%   xyrect=[xmin,xmax,ymin,ymax]: the boundaries of the stability chart
%                    (x,y):real coordinates
%   Nbasemesh=[Nx,Ny]: the size of the base mesh
%   par: system parameters
%       par.num.print: obligatory
%   stopcond(triang,points,xyrect,par): gives back the whether the just calculated triangle has a good size
%                   (0: stop, 1: separate)
%   <ctransf_x(x,direction)>: function which do the transformation in both
%                           direction, (0: x --> X & 1: X --> x, x:real ccordinate,X:inner coordinate)
%   <ctransf_y(y,direction)>: function which do the transformation in both
%                           direction, (0: y --> Y & 1: Y --> y, y:real ccordinate,Y:inner coordinate)
%   <oldres>: result from the previous calculations with zd_trianglestab to increase the resolution
%
% OUTPUT:
%   restriag.triangs
%
%   triangs: the calculated triangles
%        ID: index of the current triangle
%        pID: index of the parent triangle
%        child: indexes of the current triangle's children
%        level: its level in the family tree
%        neighbs: indexes of its point neighbourhoods
%        lineneighbs: indexes of the real neighbourhoods (only the same level)
%        points: indexes of the corner points
%        type: 0: stable, 1: unstable, 2: border
%
%   restriag.points
%
%   points: the calculated points
%        x: x coordinate of the point
%        y: y coordinate of the point
%        stab: 0: stable, 1: unstable
%        res: result of the calculation
%        triangles: indexes of the connected triangles
%        ID: index of the current point
%
%   restriag.borders
%
%   borders: the indexes of the border triangles
% VERSION: zd_trianglestab 5.5
%   1.1: pointmesh, index matrix to recognize the already calculated points.
%   1.2: some corrected error.
%   4.0: rearranged
%           23.10.2007. Corrected normalization error in other_side function
%   4.1: corrected
%   5.0: it is capable to continoue the old solution
%   5.1: par inclouded in restriag
%   5.2: bug with the ctransf_x,ctransf_y is fixed to use oldres without coordinate transformation
%   5.3: fixed to be compatible with matlab 2017
%   5.4: stoplevel has removed, stopcond giving par
%   5.5: par.num.print corrected

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATE: Vancouver, BC, Canada, UBC, 2017. 10. 09.
% Zoltan Dombovari
% Budapest University of Technology and Echonomics
% Department of Applied Mechanics
% 2007
% ddombo@gmail.com
% 2017, Ideko-IK4

%% MAIN CODE
%definition of the plotting oportunity
%  tic;
% global Ts;
% global nex;

%check whether coordinate transformations, ctransf_x
if exist('ctransf_x','var')
    if ~isempty(ctransf_x)
        exist_ctransf_x=true;
    else
        exist_ctransf_x=false;
    end
else
    exist_ctransf_x=false;
end
%check whether coordinate transformations, ctransf_y
if exist('ctransf_y','var')
    if ~isempty(ctransf_y)
        exist_ctransf_y=true;
    else
        exist_ctransf_y=false;
    end
else
    exist_ctransf_y=false;
end
%continuation of an old solution
olds=exist('oldres','var');

plotprocess=(par.num.print(2)~=0); %true: plot, false: without plot

%pointmesh to follow easily the points in its array points
%the size of the pointmesh can be determined after the first convergated triangle
f_conv_oc=false;
pID=1;
errorc=0;

TOL=1e-6;
pointmesh=[];

if ~olds


    pmdx=0;
    pmdy=0;



    if plotprocess
        %figure definition
        figure(1001);clf;
        hold on;
        dxx=(xyrect(2)-xyrect(1))*0.1;
        dyy=(xyrect(4)-xyrect(3))*0.1;
        axis([xyrect(1)-dxx xyrect(2)+dxx xyrect(3)-dyy xyrect(4)+dyy]);
        'Please arrange the windows to see the process!'
        pause;
    end

    % definition of the base datas
    xmin=xyrect(1);
    xmax=xyrect(2);
    ymin=xyrect(3);
    ymax=xyrect(4);

    Nx=Nbasemesh(1);
    Ny=Nbasemesh(2);

    dx=(xmax-xmin)/Nx;
    dy=(ymax-ymin)/Ny;

    x=xmin;
    y=ymin;
    % base mesh
    % tic;
    if errorc==0
        for iT=1:Nx+1
            for jT=1:Ny+1
                if ~exist('points','var')
                    %create the first point
                    pID=1;
                    points=cr_point(x,y,pID);
                    points((Nx+1)*(Ny+1))=points;
                else
                    pID=pID+1;
                    points(pID)=cr_point(x,y,pID);
                end
                if exist('par','var') && ~isempty(par)
                    [points(pID).res,points(pID).stab]=condfunc(trx2X(x),try2Y(y),par);
                else
                    [points(pID).res,points(pID).stab]=condfunc(trx2X(x),try2Y(y));
                end

                y=y+dy;
            end
            y=ymin;
            x=x+dx;
        end
    end
    % toc
    %create base triangles
    borders=[];
    TID=0;
    if errorc==0
        for iR=1:Nx
            for jR=1:Ny
                if ~exist('triangs','var')
                    %first triangle
                    [~,points,triangs]=cr_triang(0,1,points,[],[(iR-1)*(Ny+1)+jR+1,(iR-1)*(Ny+1)+jR,iR*(Ny+1)+jR]);
                    TID=1;
                    if triangs(TID).type==2
                        borders=[borders TID];
                    end
                    triangs(2*Nx*Ny)=triangs;
                else
                    TID=TID+1;
                    [triangs,points,triangs(TID)]=cr_triang(0,TID,points,triangs,[(iR-1)*(Ny+1)+jR+1,(iR-1)*(Ny+1)+jR,iR*(Ny+1)+jR]);
                    if triangs(TID).type==2
                        borders=[borders TID];
                    end
                end
                TID=TID+1;
                [triangs,points,triangs(TID)]=cr_triang(0,TID,points,triangs,[(iR-1)*(Ny+1)+jR+1,iR*(Ny+1)+jR+1,iR*(Ny+1)+jR]);
                if triangs(TID).type==2
                    borders=[borders TID];
                end
            end
        end
    end
else
    points=oldres.points;
    triangs=oldres.triangs;
    borders=oldres.borders;
    if plotprocess
        figure(1001);clf;
        
        zd_triangleplot(oldres.points,oldres.triangs);
    end
end
%decrease the border triangles
bid=1;

% start the fractal triangle method along the borders

zd_triangcut(bid);
% toc
% tic;
% follow the boundaries
bordID=1;
nlines=1;
IDstarts=[];
IDends=[];
direction=1; %1=+, 0=-
if errorc==0
    while bordID<length(borders) && errorc==0
        if direction==1
            % definition of a line details
            lines(nlines).IDstart=bordID;
            IDstarts=[IDstarts borders(bordID)];
            new_line_tr=1;
        end
        km=1;
        wrong_TID=12;
        while wrong_TID~=0 && errorc==0
            found=10;
            newline=false;
            while found~=0
                %just something big number nothing else
                mpoints=2^16;
                if plotprocess
                    plot_triang(points,triangs(borders(bordID)),'k');
                end
                %neighborhood location

                found=0;

                for mN=1:length(triangs(borders(bordID)).lineneighbs)
                    if bordID>1
                        %find the neighbourhood in the borders
                        old_neighbours=find(borders(1:bordID-1)==triangs(borders(bordID)).lineneighbs(mN));
                    else
                        old_neighbours=[];
                    end
                    %if it is exist, it must be the same type and level
                    if (isempty(old_neighbours))&&...
                            triangs(triangs(borders(bordID)).lineneighbs(mN)).level==triangs(borders(bordID)).level &&...
                            triangs(triangs(borders(bordID)).lineneighbs(mN)).type==triangs(borders(bordID)).type
                        found=mN;
                        %investigation of the line boundaries (start)
                        if ~isempty(IDstarts)
                            id=find(IDstarts==triangs(triangs(borders(bordID)).lineneighbs(found)).ID);
                            if ~isempty(id)
                                found=0;
                                newline=true;
                            end
                        end
                        %(end)
                        if ~isempty(IDends) && found~=0
                            id=find(IDends==triangs(triangs(borders(bordID)).lineneighbs(found)).ID);
                            if ~isempty(id)
                                found=0;
                                newline=true;
                            end
                        end
                    end
                end

                if found~=0
                    id3=triangs(borders(bordID)).lineneighbs(found);
                    if plotprocess
                        plot_triang(points,triangs(triangs(borders(bordID)).lineneighbs(found)),'y');
                        pause(0.01);
                    end
                    id=find(borders==id3);
                    if ~isempty(id)
                        %define the next step
                        bordID=bordID+1;
                        if plotprocess && bordID>length(borders)
                            for iP=1:length(borders)
                                plot_triang(points,triangs(borders(iP)),'m');
                            end
                        end
                        new_line_tr=new_line_tr+1;
                        id2=borders(id);
                        borders(id)=borders(bordID);
                        borders(bordID)=id2;
                    else
                        error('missing triangle');
                    end
                end
            end

            km=km+1;
            if plotprocess
                plot_triang(points,triangs(borders(bordID)),'m');
            end
            wrong_TID=0;
            %if is not exist a real neighbourhood
            if found==0 && ~newline
                if length(triangs(borders(bordID)).lineneighbs)<3
                    pt=[0 0 0]';
                    for iT=1:length(triangs(borders(bordID)).lineneighbs)
                        for jT=1:3
                            id=find(triangs(borders(bordID)).points==triangs(triangs(borders(bordID)).lineneighbs(iT)).points(jT));
                            if ~isempty(id)
                                pt(id)=pt(id)+1;
                            end
                        end
                    end
                    % find the side of the triangle where it has not got real
                    % any type of neighbours in the same level
                    [~,mis_neigh_pt]=min(pt);
                    [~,mis_neigh_pt(2)]=min([pt(1:mis_neigh_pt(1)-1); 1000; pt(mis_neigh_pt(1)+1:end)]);
                    akt_triangID=borders(bordID);
                    found=10;
                    % find the closest non real neighbourhood in the lower
                    % levels
                    while found~=0 && triangs(akt_triangID).level~=0

                        if found~=0 && length(triangs(triangs(akt_triangID).pID).lineneighbs)>=1
                            for kB=1:length(triangs(triangs(akt_triangID).pID).lineneighbs)
                                res=other_side(points,triangs(borders(bordID)),mis_neigh_pt,triangs(triangs(triangs(akt_triangID).pID).lineneighbs(kB)),xyrect);
                                if res==1
                                    found=0;
                                    wrong_TID=triangs(triangs(akt_triangID).pID).lineneighbs(kB);
                                end
                            end

                            akt_triangID=triangs(akt_triangID).pID;
                        end
                    end
                    if wrong_TID~=0
                        if plotprocess
                            plot_triang(points,triangs(wrong_TID),'c');
                        end
                        borders=[borders wrong_TID];
                        id4=length(borders);
                        % cut the non real neighbourhood
                        zd_triangcut(length(borders));
                        if id4==length(borders) && triangs(borders(end)).type~=2
                            borders=borders(1:end-1);
                        else
                            bordID=bordID-1;
                        end
                    end
                else
                    if plotprocess
                        plot_triang(points,triangs(borders(bordID)),'y');
                    end
                end
            end
            if wrong_TID==0
                %if there is not any real and non real neighbourhoods, then it
                %is probably the border of the stability chart, change the
                %direction
                if direction==0
                    lines(nlines).IDend=bordID;
                    IDends=[IDends borders(bordID)];
                    bordID=bordID+1;
                    nlines=nlines+1;
                    direction=1;
                else
                    direction=0;
                    borders(lines(nlines).IDstart:bordID)=borders(bordID:-1:lines(nlines).IDstart);
                end
            end
        end
    end
end
% toc

restriag.triangs=triangs;
restriag.points=points;
restriag.borders=borders;

restriag.par=par;

return;

    function zd_triangcut(bid)

        % function [points,triangs,borders]=zd_triangcut(condfunc,points,triangs,borders,bid,stopcond)
        %  "the representation of the triangle fractals method"
        % INPUT:
        %   condfunc,xyrect,Nbasemesh,stopcond: description in the heading
        %   points: already calculated points
        %   triangs: already calculated triangles
        %   bid: position in the borders dynamic array
        % OUTPUT:
        %   triangs,points,borders: definition in the heading

        % Zoltan Dombovari
        % Budapest University of Technology and Echonomics
        % Department of Applied Mechanics
        % 2007
        % ddombo@gmail.com
        % tic;
        pID=length(points);
        TID=length(triangs);
        % scan the border triangles
        while bid<=length(borders)
            % investigation of the stop condition (external function)
            if stopcond(triangs(borders(bid)),points,xyrect,par)==1
                %definition of the new point
                Npx=(points(triangs(borders(bid)).points(1)).x+points(triangs(borders(bid)).points(3)).x)/2;
                Npy=(points(triangs(borders(bid)).points(1)).y+points(triangs(borders(bid)).points(3)).y)/2;
                %investigation is it an existing points?
                res=p_exist(points,Npx,Npy,TOL);
                if res==0
                    %new point if not exist
                    pID=pID+1;
                    points(pID)=cr_point(Npx,Npy,pID);
                    if exist('par','var') && ~isempty(par)
                        [points(pID).res,points(pID).stab]=condfunc(trx2X(Npx),try2Y(Npy),par);
                    else
                        [points(pID).res,points(pID).stab]=condfunc(trx2X(Npx),try2Y(Npy));
                    end
                    CpID=pID;
                else
                    %index of the existed point if exist
                    CpID=res;
                end
                %the new children matrixes
                TID=TID+1;
                [triangs,points,triangs(TID)]=cr_triang(triangs(borders(bid)).ID,TID,points,triangs,[triangs(borders(bid)).points(1),CpID,triangs(borders(bid)).points(2)]);
                TID=TID+1;
                [triangs,points,triangs(TID)]=cr_triang(triangs(borders(bid)).ID,TID,points,triangs,[triangs(borders(bid)).points(3),CpID,triangs(borders(bid)).points(2)]);
                %if the first triangle is border, it is put the place of its parent

                if triangs(TID-1).type==2
                    borders(bid)=TID-1;
                end
                if triangs(TID).type==2
                    if borders(bid)==TID-1
                        %if the second also border it is put to an extended place
                        borders=[borders TID];
                    else
                        %or not
                        borders(bid)=TID;
                    end
                elseif triangs(TID-1).type~=2
                    %if neither one are border type step to the next border triangle
                    bid=bid+1;
                end
            else
                % if the stopcond is granted take the next border triangle
                if ~f_conv_oc
                    pointmesh=create_pointmesh(points,triangs,xyrect,triangs(borders(bid)).ID);
                    f_conv_oc=true;
                end
                bid=bid+1;
            end
        end
        % Ts=Ts+toc;
    end

    function point=cr_point(x,y,pID)
        % function point=cr_point(x,y)
        %  "create a point structure"
        % INPUT:
        %   x: x coordinate of the point
        %   y: y coordinate of the point
        % OUTPUT:
        %   point: created point structure, see the heading
        %   pointmesh:???

        % Zoltan Dombovari
        % Budapest University of Technology and Echonomics
        % Department of Applied Mechanics
        % 2007
        % ddombo@gmail.com
        % tic;

        point.x=x;
        point.y=y;
        point.stab=0;
        point.trlevstore=-1;
        point.triangles=[];
        point.res=[];
        point.ID=pID;
        if ~isempty(pointmesh)
            pointmesh=put_ptopmesh(point,pointmesh);
        end
        % Ts=Ts+toc;
    end

    function [triangs,points,triang]=cr_triang(parentID,TID,points,triangs,indexes);
        % function [triangs,points,triang]=cr_triang(parentID,TID,points,triangs,indexes);
        %  "create a new triangle structure"
        % INPUT:
        %   parentID: index of parent triangle of the new triangle (0: base mesh:  triangle, ancestor triangle)
        %   TID: index of the new triangle
        %   points, triangs: see the heading
        %   indexes: of the corner points of the new triangle
        % OUTPUT:
        %   points, triangs: see the heading
        %   triang: structure of the new triang
        % Zoltan Dombovari
        % Budapest University of Technology and Echonomics
        % Department of Applied Mechanics
        % 2007
        % ddombo@gmail.com
        triang.ID=TID;
        triang.pID=parentID;
        %         triang.child=[0 0];
        if parentID~=0
            triang.level=triangs(parentID).level+1;
            %             if triangs(parentID).child(1)==0
            %                 triangs(parentID).child(1)=TID;
            %             else
            %                 triangs(parentID).child(2)=TID;
            %             end
        else
            triang.level=0;
        end

        triang.neighbs=[];
        triang.lineneighbs=[];
        triang.points=indexes;
        %neighbourhood investigation
        %          tic;
        for k=1:3
            if ~isempty(triangs)
                if ~isempty(points(indexes(k)).triangles) && points(indexes(k)).trlevstore(1)>=triang.level

                    for l=1:points(indexes(k)).trlevstore(triang.level+2)
                        idtr=points(indexes(k)).triangles(triang.level+1,l);
                        %                         [idtr length(triangs)]

                        if triangs(idtr).ID~=triang.ID

                            if isempty(find(triang.neighbs==triangs(points(indexes(k)).triangles(triang.level+1,l)).ID));
                                triang.neighbs=[triang.neighbs triangs(points(indexes(k)).triangles(triang.level+1,l)).ID];
                                triangs(points(indexes(k)).triangles(triang.level+1,l)).neighbs=[triangs(points(indexes(k)).triangles(triang.level+1,l)).neighbs triang.ID];
                            else
                                id_com_points=0;
                                for m=1:3
                                    id2=find(triang.points==triangs(points(indexes(k)).triangles(triang.level+1,l)).points(m));
                                    if ~isempty(id2)
                                        id_com_points=id_com_points+1;
                                    end
                                end
                                if id_com_points==2 && isempty(find(triang.lineneighbs==triangs(points(indexes(k)).triangles(triang.level+1,l)).ID))
                                    triang.lineneighbs=[triang.lineneighbs triangs(points(indexes(k)).triangles(triang.level+1,l)).ID];
                                    triangs(points(indexes(k)).triangles(triang.level+1,l)).lineneighbs=[triangs(points(indexes(k)).triangles(triang.level+1,l)).lineneighbs triang.ID];
                                end
                            end
                        end
                    end
                end
            end
        end
        for k=1:3
            %put the new triangles in the corner points triangle list
            %the structure of triangles
            %pr lev store [No tr.level | No stored tr. on lev0 | No stored tr. on lev1|...|No stored tr. on levp]
            %level 0 [tr1         | tr2                   |...]
            %level 1 [tr1      ...
            if points(indexes(k)).trlevstore(1)>=triang.level
                points(indexes(k)).trlevstore(triang.level+2)=points(indexes(k)).trlevstore(triang.level+2)+1;
            else
                points(indexes(k)).trlevstore(1)=triang.level;
                points(indexes(k)).trlevstore(triang.level+2)=1;
            end
            points(indexes(k)).triangles(triang.level+1,points(indexes(k)).trlevstore(triang.level+2))=TID;
        end
        %          toc;
        %estabilish the type of the new triangle
        if points(indexes(1)).stab==points(indexes(2)).stab && points(indexes(2)).stab==points(indexes(3)).stab
            triang.type=points(indexes(1)).stab;
        else
            triang.type=2;
        end

        if plotprocess
            switch triang.type
                case {0}
                    style='g';
                case {1}
                    style='r';
                case {2}
                    style='b';
            end
            plot_triang(points,triang,style);
            %             sx=0;
            %             sy=0;
            %             for l=1:3
            %                 sx=sx+points(triang.points(l)).x/3;
            %                 sy=sy+points(triang.points(l)).y/3;
            %             end
            %             text(sx,sy,['---' int2str(triang.ID)]);
            %             if triang.ID==58
            %                 plot(points(triang.points(1)).x,points(triang.points(1)).y,'mo');
            %                 plot(points(triang.points(2)).x,points(triang.points(2)).y,'mv');
            %                 plot(points(triang.points(3)).x,points(triang.points(3)).y,'m*');
            %             end
            pause(0.01);

        end
    end

    function res=p_exist(points,Nx,Ny,TOL)
        % function res=p_exist(points,Nx,Ny,TOL)
        %  "investigation of the points, is it exist already?"
        % INPUT:
        %   points: see the heading
        %   Nx,Ny: supposed points coordinates
        %   TOL: tolerance of the equality
        % OUTPUT:
        %   res: 0: not the same,  index: if the same
        % Zoltan Dombovari
        % Budapest University of Technology and Echonomics
        % Department of Applied Mechanics
        % 2007
        % ddombo@gmail.com

        %          tic;
        lx=xyrect(2)-xyrect(1);
        ly=xyrect(4)-xyrect(3);
        if ~isempty(pointmesh)
            a=ceil((Nx-xyrect(1)+TOL)/pmdx);
            b=ceil((Ny-xyrect(3)+TOL)/pmdy);
            if ~isempty(pointmesh(a,b).points)
                for j=1:length(pointmesh(a,b).points)

                    id2(j)=(((points(pointmesh(a,b).points(j)).x-Nx)/lx).^2+((points(pointmesh(a,b).points(j)).y-Ny)/ly).^2);
                end
                id=find(id2<TOL);
                if isempty(id)
                    res=0;
                else
                    res=pointmesh(a,b).points(id);
                    if plotprocess
                        % plot(Nx,Ny,'mo');
                    end
                end
            else
                res=0;
            end
        else
            res=0;
        end
        %          Ts=Ts+toc;
    end

    function res=other_side(points,ortriang,pt2,newtriang,xyrect);
        % function res=other_side(points,ortriang,pt2,newtriang);
        %  "investigate the non real neighborhoods "
        % INPUT:
        %   points: see the heading
        %   ortriang: structure of the original triangle
        %   pt2: its two corner points which define a line
        %   newtriang: investigated triangle
        % OUTPUT:
        %   res: 1: non real neighbourhood, 0: non neighbourhood
        % Zoltan Dombovari
        % Budapest University of Technology and Echonomics
        % Department of Applied Mechanics
        % 2007
        % ddombo@gmail.com
        % tic;
        lx=xyrect(2)-xyrect(1);
        ly=xyrect(4)-xyrect(3);

        ors1.x=0;
        ors1.y=0;
        news2.x=0;
        news2.y=0;
        for i=1:3
            ors1.x=ors1.x+(points(ortriang.points(i)).x-xyrect(1))/lx;
            ors1.y=ors1.y+(points(ortriang.points(i)).y-xyrect(3))/ly;
            news2.x=news2.x+(points(newtriang.points(i)).x-xyrect(1))/lx;
            news2.y=news2.y+(points(newtriang.points(i)).y-xyrect(3))/ly;
        end

        ors1.x=ors1.x/3;
        ors1.y=ors1.y/3;

        news2.x=news2.x/3;
        news2.y=news2.y/3;

        v1.x=((points(ortriang.points(pt2(2))).x-xyrect(1))-(points(ortriang.points(pt2(1))).x-xyrect(1)))/lx;
        v1.y=((points(ortriang.points(pt2(2))).y-xyrect(3))-(points(ortriang.points(pt2(1))).y-xyrect(3)))/ly;

        v2.x=ors1.x-news2.x;
        v2.y=ors1.y-news2.y;

        if v2.x*v1.y-v2.y*v1.x~=0
            p2y=(points(ortriang.points(pt2(1))).y-xyrect(3))/ly;
            p2x=(points(ortriang.points(pt2(1))).x-xyrect(1))/lx;
            v=(v1.x*(news2.y-p2y)-v1.y*(news2.x-p2x))/(v2.x*v1.y-v2.y*v1.x);
            if v>0 && v<1
                res=1;
            else
                res=0;
            end
        else
            res=0;
            errorc=1;
            'errorc'
        end

        % Ts=Ts+toc;
    end

    function pointmesh=create_pointmesh(points,triangs,xyrect,f_convtrID)
        % tic;
        for i=1:3
            xp(i)=points(triangs(f_convtrID).points(i)).x;
            yp(i)=points(triangs(f_convtrID).points(i)).y;
        end
        pmdx=(max(xp)-min(xp))*2;
        pmdy=(max(yp)-min(yp))*2;
        Nx=ceil((xyrect(2)-xyrect(1)+TOL)/pmdx);
        Ny=ceil((xyrect(4)-xyrect(3)+TOL)/pmdy);
        pointmesh(1).points=[];
        pointmesh(Nx,Ny)=pointmesh(1);
        for i=2:length(points)
            pointmesh=put_ptopmesh(points(i),pointmesh);
        end
        %         Ts=Ts+toc;
    end

    function pointmesh=put_ptopmesh(point,pointmesh)
        %  tic;
        a=ceil((point.x-xyrect(1)+TOL)/pmdx);
        b=ceil((point.y-xyrect(3)+TOL)/pmdy);
        pointmesh(a,b).points=[pointmesh(a,b).points point.ID];
        %         Ts=Ts+toc;
    end

    function resX=trx2X(x)
        if exist_ctransf_x
            resX=ctransf_x(x,0);
        else
            resX=x;
        end
    end

    function resx=trX2x(X)
        if exist_ctransf_x
            resx=ctransf_x(x,1);
        else
            resx=X;
        end
    end

    function resY=try2Y(y)
        if exist_ctransf_y
            resY=ctransf_y(y,0);
        else
            resY=y;
        end
    end

    function resy=trY2y(y)
        if exist_ctransf_y
            resY=ctransf_y(y,0);
        else
            resY=y;
        end
    end

end