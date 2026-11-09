function plot_stab_mu(restriang,p_type)
%PLOT_STAB_lambda Plot stability map based on output of fractal search 
% -- assuming res.aceig stores critical multipliers --
%   restriang: output data structure from fractal method
%   p_type: output plot: 0) line, 1) contour, 2) heatmap, 3) frequency, 4) scatter 

if nargin<2
    p_type = 1; % create simple contour plot by default
end

%select the points
X = zeros(1,length(restriang.points)); Y = X; Z = X; f = X;
for i=1:length(restriang.points)
    X(i) = restriang.points(i).x;
    Y(i) = restriang.points(i).y;
    Z(i) = restriang.points(i).res.aceig;
    if p_type == 3
        f(i) = restriang.points(i).res.om_c;
    end
end

% grid with minimal required resolution
dx_min = min(diff(unique(X,'sorted')));
dy_min = min(diff(unique(Y,'sorted')));
xi = min(X):dx_min/2:max(X);
yi = min(Y):dy_min/2:max(Y);
[Xi,Yi] = meshgrid(xi,yi);

%interpolation of the stability
xy = (max(X)-min(X))/(max(Y)-min(Y));
Zi = griddata(X,xy*Y,Z,Xi,xy*Yi,'linear');

switch p_type
    case 0
        contour(Xi,Yi,Zi,[1 1],'k','LineWidth',1.5);
    case 1
        % contour plot
        contourf(Xi,Yi,Zi,[0 1 1]); hold on
        contour(Xi,Yi,Zi,[1 1],'LineWidth',2,'EdgeColor','k')
        colormap([0.9 0.9 0.9; 1.0 1.0 1.0])
    case 2
        % heatmap of characteristic multipliers
        contourf(Xi,Yi,Zi,50,'EdgeColor','none'); hold on
        contour(Xi,Yi,Zi,[1 1],'LineWidth',2,'EdgeColor','w')
        c = colorbar;
        c.Label.Interpreter = 'latex';
        c.Label.String = '$\Re\{\lambda_{\mathrm{crit}}\}$';
        hold off
    case 3
        % dominant chatter frequencies
        c = contour(Xi,Yi,Zi,[1 1]);
        fi = griddata(X,Y,f,Xi,Yi,'linear');
        fc = zeros(1,size(c,2));
        for i=1:size(c,2)
            if Xi(1,1)<=c(1,i) && c(1,i)<=Xi(1,end) && Yi(1,1)<=c(2,i) && c(2,i)<=Yi(end,1)
                fc(i)=interp2(Xi,Yi,fi,c(1,i),c(2,i));
            else
                fc(i)=0;
            end
        end
        fc(c(1,:)>max(xi)) = NaN;
        fc(c(1,:)<min(xi)) = NaN;
        plot(c(1,2:end),fc(2:end),'k','Linewidth',2);
    case 4
        % scatter map for validation
        scatter(X(Z>1),Y(Z>1),'rx'); hold on
        scatter(X(Z<=1),Y(Z<=1),'go')
        ylim([min(Y),max(Y)]);
        box on
end
xlim([min(X),max(X)]);
hold off
end

