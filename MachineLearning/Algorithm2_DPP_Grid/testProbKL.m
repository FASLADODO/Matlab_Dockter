nn = 100;

p = linspace(0.0001,1,nn)';

[X,Y] = meshgrid(p,p);
Z = log10(X./Y);

figure
mesh(X,Y,Z);