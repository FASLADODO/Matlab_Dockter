
nn = 100;

p = linspace(0.001,1,nn)';

[X,Y] = meshgrid(p,p);

%functions
Z1 = log10(X./Y);
Z2 = exp(-X./Y);
Z3 = X.*log10(X./Y);
Z4 = X.*log10(X./Y) + Y.*log10(Y./X);
Z5 = (X+Y).*log(X./Y);

% figure
% mesh(X,Y,Z1);
% title('log')
% 
% figure
% mesh(X,Y,Z2);
% title('exp')
% 
% figure
% mesh(X,Y,Z3);
% title('p log p')
% 
% figure
% mesh(X,Y,Z4);
% title('p log p + p log p ')

figure
mesh(X,Y,Z5);
%title('p+p log p ')
xlabel('P(C_{1})','FontSize',fsize)
ylabel('P(C_{2})','FontSize',fsize)
zlabel('W_{KL}','FontSize',fsize)
colormap(flipud(cool))
hc = colorbar;
ylabel(hc, 'W_{KL}','FontSize',fsize)

%% for combining KL and density\

fsize = 14;

nn = 101;

w = linspace(-3,3,nn)';
w(w==0) = [];
p = linspace(0,1,nn)';
p(p==0) = [];

[X,Y] = meshgrid(p,w);

Z1 = X .* Y;
Z2 = ((X .* Y) ./ (abs(X) + abs(Y)));
Z3 = X .* norm(log10(Y));

figure
mesh(X,Y,Z1);
title('product')

figure
mesh(X,Y,Z2);
xlabel('\rho','FontSize',fsize)
ylabel('W_{KL}','FontSize',fsize)
zlabel('S_{j,k}','FontSize',fsize)
colormap(flipud(cool))
hc = colorbar;
ylabel(hc, 'S_{j,k}','FontSize',fsize)


%% other way

tt1 = 0:0.001:1;
tt2 = 1-tt1;

figure
scatter3(tt1,tt2,(tt1+tt2).*log10(tt1./tt2))