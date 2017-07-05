function handle = Surface3Dalt(x,y,z,limits,nn)
%Surface3D Summary: plot 3D surface from x,y,z vectors

%linear spacing
xlin = linspace(limits(1),limits(2),nn);
ylin = linspace(limits(3),limits(4),nn);

%mesh grid
[X,Y] = meshgrid(xlin,ylin);
f = scatteredInterpolant(x,y,z);
Z = f(X,Y);

%plot
handle = mesh(X,Y,Z);

end

