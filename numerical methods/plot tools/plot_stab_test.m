function plot_stab_test(restriang,p_type)
%PLOT_STAB_lambda Plot stability map based on output of fractal search
% -- assuming res.aceig stores charachteristic exponents --
%   restriang: output data structure from fractal method
%   p_type: output plot: 0) stab 1) aceig, 2) feig
if nargin<2
    p_type = 0; % plot stability by default
end

%select the points
X = zeros(1,length(restriang.points)); Y = X; Z = X; f = X;
for i=1:length(restriang.points)
    X(i) = restriang.points(i).x;
    Y(i) = restriang.points(i).y;
    switch p_type
        case 0
            Z(i) = restriang.points(i).stab;
        case 1
            Z(i) = sign(restriang.points(i).res.aceig);
        case 2
            Z(i) = sign(restriang.points(i).res.feig);
    end
end

% Scatter plot
scatter(X(Z>0),Y(Z>0),'rx'); hold on; box on
scatter(X(Z<=0),Y(Z<=0),'go'); hold off
axis([min(X),max(X),min(Y),max(Y)]);

end

