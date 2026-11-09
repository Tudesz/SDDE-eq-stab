function plot_triang(points,triang,style);

for l=1:3
    kx(l)=points(triang.points(l)).x;
    ky(l)=points(triang.points(l)).y;
end;
fill(kx,ky,style);