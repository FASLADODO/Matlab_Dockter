load fisheriris

columns = [3,4]; %PL and PW

X = meas(:,columns);

%build classifier
MdlLinear = fitcdiscr(X,species);

%classes
classes = unique(species)

d1 = meas( strcmp(species,classes(1)) ,columns);
d2 = meas( strcmp(species,classes(2)) ,columns);
d3 = meas( strcmp(species,classes(3)) ,columns);

%coefficients from classifier
K12 = MdlLinear.Coeffs(1,2).Const;
L12 = MdlLinear.Coeffs(1,2).Linear;
K23 = MdlLinear.Coeffs(2,3).Const;
L23 = MdlLinear.Coeffs(2,3).Linear;


figure
h1 = scatter(d1(:,1),d1(:,2),'r+')
hold on
h2 = scatter(d2(:,1),d2(:,2),'b*')
hold on
h3 = scatter(d3(:,1),d3(:,2),'gd')

hold on

%plot boundary
nn =25;
xp = [ linspace(1,7,nn)' ];
f12 = -(xp*L12(1) + repmat(K12,nn,1))./L12(2);
f23 = -(xp*L23(1) + repmat(K23,nn,1))./L23(2);

h4 = plot(xp,f12,'b-')
hold on
h5 = plot(xp,f23,'g-')
hold off

legend('Setosa','Versicolor','Virginica','Location','best')

%% linear 3D
NN = 100;

%create some rotated gaussian data
cv1 = [2,0.1,0.1;
    0.1,2.5,0.1;
    0.1,0.1,3];
mu1 = [3,4,5];

cv2 = [2,0.1,0.1;
    0.1,2.5,0.1;
    0.1,0.1,3];
mu2 = [1,-1,1];


d1 = mvnrnd(mu1,cv1,NN);
d2 = mvnrnd(mu2,cv2,NN);
D = [d1;d2];
labels = [ones(NN,1);ones(NN,1)*2];

MdlLin3D = fitcdiscr(D,labels);


figure
scatter3(d1(:,1),d1(:,2),d1(:,3),'r+')
hold on
scatter3(d2(:,1),d2(:,2),d2(:,3),'b*')
hold on
PlotLDA3D(MdlLin3D,D,25);
hold off

%% quadratic 3d
NN = 100;

%create some rotated gaussian data
cv1 = [2,0.1,0.1;
    0.1,2.5,0.1;
    0.1,0.1,3];
mu1 = [3,4,5];

cv2 = [2,0.1,0.1;
    0.1,2.5,0.1;
    0.1,0.1,3];
mu2 = [1,-1,1];


d1 = mvnrnd(mu1,cv1,NN);
d2 = mvnrnd(mu2,cv2,NN);
D = [d1;d2];
labels = [ones(NN,1);ones(NN,1)*2];

MdlQuadratic = fitcdiscr(D,labels,'DiscrimType','quadratic');

nn = 25;

%     lims = [min(D); max(D)];
% 
%     K = MdlQuadratic.Coeffs(1,2).Const;
%     L = MdlQuadratic.Coeffs(1,2).Linear;
%     Q = MdlQuadratic.Coeffs(1,2).Quadratic;
% 
%     %space out some vectors
%     xv = linspace(lims(1,1),lims(2,1),nn);
%     yv = linspace(lims(1,2),lims(2,2),nn);
%     zv = linspace(lims(1,3),lims(2,3),nn);
%     [xx,yy,zz] = meshgrid(xv,yv,zv);
%     
%     %functions
%     f = @(x,y,z) K + [x y z]*L + sum([x y z] .* ([x y z]*Q), 2);
%     %reshapce
%     v = f(xx(:),yy(:),zz(:));
%     v = reshape(v,size(xx));
    
figure
scatter3(d1(:,1),d1(:,2),d1(:,3),'r+')
hold on
scatter3(d2(:,1),d2(:,2),d2(:,3),'b*')
hold on
PlotQDA3D(MdlQuadratic,D,25);
hold off

[plabel score cost] = predict(MdlQuadratic,D);
