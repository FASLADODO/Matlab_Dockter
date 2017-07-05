%make a huge pile of data

dist = 3;
nn = [500,500];

SS = 8; %dimensions

means{1} = ones(1,SS);
sigmas{1} = eye(SS);

means{2} = means{1} + ones(1,SS)*dist/sqrt(SS);
sigmas{2} = eye(SS);

[Data,Labels] = CreateData(means,sigmas,nn);

testcolumns = 1:8;
columnlabels = {'1','2','3','4','5','6','7','8'};

[Variations,BestSep,BestSepClass,Data_All,Diff_All] = RelativeRBFSeperabilityCheck(Data,Labels,testcolumns,columnlabels,2,false);
    

BestSep.bestnumcol
BestSep.bestcolumns

%% Make some 1D rand data

dist = 3;
SS = 1;

means{1} = ones(1,SS)*2;
sigmas{1} = eye(SS)./sqrt(SS); %it still scales down with this

means{2} = means{1} + ones(1,SS).*dist/sqrt(SS);
sigmas{2} = eye(SS)./sqrt(SS);

nn = [5000,5000];

[Data,Labels] = CreateData(means,sigmas,nn);

sig = norm(std(Data));
n = length(Data);
bw = 1.06*sig*(n^(-1/5))

figure
gscatter(Data(:,1),zeros(length(Data),1),Labels,'br')

%compute dem rbfs
D1 = Data(Labels == 1,:);
D2 = Data(Labels == 2,:);
[pwithin1,pbetween1] = ComputeDiscriminateRBF(D1,D2,bw);
[pwithin2,pbetween2] = ComputeDiscriminateRBF(D2,D1,bw);

BD1 = BhattacharyyaDistance(pwithin1,pbetween1)
BD2 = BhattacharyyaDistance(pwithin2,pbetween2)

%check scaling
[scale1,Grid1,p1] = integralRBF(D1,bw);
[scale2,Grid2,p2] = integralRBF(D2,bw);
% 
% scalew1 = ScaleRBF(D1,pwithin1)
% scalew2 = ScaleRBF(D2,pwithin2)
% 
% Q1 = trapz(Grid1,p1)
% Q2 = trapz(Grid2,p2)


figure
gscatter(Data(:,1),ones(length(Data),1)*-0.001,Labels,'br')
hold on
scatter(Grid1(:,1),p1,'g.'); 
hold on
scatter(Grid2(:,1),p2,'g.'); 
hold off


Difference1 = ComputeRBFDifference(pwithin1,pbetween1);
Difference2 = ComputeRBFDifference(pwithin2,pbetween2);


figure
gscatter(Data(:,1),ones(length(Data),1)*-0.001,Labels,'br')
hold on
h1 = scatter(D1(:,1),pwithin1,'g.');
hold on
h2 = scatter(D1(:,1),pbetween1,'c.');
hold off
xlabel('samples')
ylabel('probabilities')
legend([h1(1),h2(1)],'within1','between1')
title('probabilities scaled class 1')

figure
gscatter(Data(:,1),ones(length(Data),1)*-0.001,Labels,'br')
hold on
h1 = scatter(D2(:,1),pwithin2,'g.');
hold on
h2 = scatter(D2(:,1),pbetween2,'c.');
hold off
xlabel('samples')
ylabel('probabilities')
legend([h1(1),h2(1)],'within2','between2')
title('probabilities scaled class 2')

figure
gscatter(Data(:,1),ones(length(Data),1)*-0.001,Labels,'br')
hold on
h1 = scatter(D1(:,1),Difference1,'m.');
hold on
h2 = scatter(D2(:,1),Difference2,'y.');
hold off
xlabel('samples')
ylabel('Seperability')
legend([h1(1),h2(1)],'Relative1','Relative2')
title('Seperability RBF scaled')

%% test 2D data

SS = 2;

means{1} = ones(1,SS)*2;
sigmas{1} = eye(SS)./sqrt(SS);

means{2} = means{1} + ones(1,SS).*dist/sqrt(SS);
sigmas{2} = eye(SS)./sqrt(SS);

[Data,Labels] = CreateData(means,sigmas,nn);

figure
gscatter(Data(:,1),Data(:,2),Labels,'br');
title('initial data')

sig = norm(std(Data));
n = length(Data);
bw = 1.06*sig*(n^(-1/5))


%compute dem rbfs
D1 = Data(Labels == 1,:);
D2 = Data(Labels == 2,:);
[pwithin1,pbetween1] = ComputeDiscriminateRBF(D1,D2,bw);
[pwithin2,pbetween2] = ComputeDiscriminateRBF(D2,D1,bw);

Difference1 = ComputeRBFDifference(pwithin1,pbetween1);
Difference2 = ComputeRBFDifference(pwithin2,pbetween2);

BD1 = BhattacharyyaDistance(pwithin1,pbetween1)
BD2 = BhattacharyyaDistance(pwithin2,pbetween2)

figure
gscatter(Data(:,1),Data(:,2),Labels,'br');
hold on
h1 = Surface3D(D1(:,1),D1(:,2),pwithin1);
hold on
h2 = Surface3D(D1(:,1),D1(:,2),pbetween1);
hold off
xlabel('dim1')
ylabel('dim2')
zlabel('probabilities')
legend([h1(1),h2(1)],'within1','between1')
title('probabilities scaled class 1')
view(45,45)

figure
gscatter(Data(:,1),Data(:,2),Labels,'br');
hold on
h1 = Surface3D(D2(:,1),D2(:,2),pwithin2);
hold on
h2 = Surface3D(D2(:,1),D2(:,2),pbetween2);
hold off
xlabel('dim1')
ylabel('dim2')
zlabel('probabilities')
legend([h1(1),h2(1)],'within2','between2')
title('probabilities scaled class 2')
view(45,45)

figure
gscatter(Data(:,1),Data(:,2),Labels,'br');
hold on
h1 = Surface3D(D1(:,1),D1(:,2),Difference1);
hold on
h2 = Surface3D(D2(:,1),D2(:,2),Difference2);
hold off
xlabel('dim1')
ylabel('dim2')
zlabel('Seperability')
legend([h1(1),h2(1)],'Relative1','Relative2')
title('Seperability RBF scaled')
view(45,45)



%% test 3D data

plotsize3 = 2;

SS = 3;

means{1} = ones(1,SS)*2;
sigmas{1} = eye(SS)./sqrt(SS);

means{2} = means{1} + ones(1,SS).*dist/sqrt(SS);
sigmas{2} = eye(SS)./sqrt(SS);

[Data,Labels] = CreateData(means,sigmas,nn);


sig = norm(std(Data));
n = length(Data);
bw = 1.06*sig*(n^(-1/5))

figure
gscatter3(Data(:,1),Data(:,2),Data(:,3),Labels,'br');
title('initial data')


%compute dem rbfs
D1 = Data(Labels == 1,:);
D2 = Data(Labels == 2,:);
[pwithin1,pbetween1] = ComputeDiscriminateRBF(D1,D2,bw);
[pwithin2,pbetween2] = ComputeDiscriminateRBF(D2,D1,bw);

Difference1 = ComputeRBFDifference(pwithin1,pbetween1);
Difference2 = ComputeRBFDifference(pwithin2,pbetween2);

BD1 = BhattacharyyaDistance(pwithin1,pbetween1)
BD2 = BhattacharyyaDistance(pwithin2,pbetween2)

figure
h1 = scatter3(D1(:,1),D1(:,2),D1(:,3),plotsize3,pwithin1);
hold on
h2 = scatter3(D2(:,1),D2(:,2),D2(:,3),plotsize3,pwithin2);
hold off
xlabel('dim1')
ylabel('dim2')
zlabel('dim3')
legend([h1(1),h2(1)],'within1','within2')
title('probabilities scaled within both classes')
colormap cool
colorbar
view(45,45)

figure
h1 = scatter3(D1(:,1),D1(:,2),D1(:,3),plotsize3,pbetween1);
hold on
h2 = scatter3(D2(:,1),D2(:,2),D2(:,3),plotsize3,pbetween2);
hold off
xlabel('dim1')
ylabel('dim2')
zlabel('dim3')
legend([h1(1),h2(1)],'between1','between2')
title('probabilities scaled between both classes')
colormap cool
colorbar
view(45,45)

figure
h1 = scatter3(D1(:,1),D1(:,2),D1(:,3),plotsize3,Difference1);
hold on
h2 = scatter3(D2(:,1),D2(:,2),D2(:,3),plotsize3,Difference2);
hold off
xlabel('dim1')
ylabel('dim2')
zlabel('dim3')
legend([h1(1),h2(1)],'Relative1','Relative2')
title('Seperability RBF scaled')
colormap cool
colorbar
view(45,45)


%% loop through dims, figure out relationship to max diff

dimensionz = 1:8;
nn = [5000,5000];

storeit = [];

for SS = dimensionz

    means{1} = ones(1,SS)*2;
    sigmas{1} = eye(SS);

    means{2} = means{1} + ones(1,SS)*dist/sqrt(SS);
    sigmas{2} = eye(SS);

    [Data,Labels] = CreateData(means,sigmas,nn);
    
    sig = norm(std(Data));
    n = length(Data);
    bw = 1.06*sig*(n^(-1/5));

    %compute dem rbfs
    D1 = Data(Labels == 1,:);
    D2 = Data(Labels == 2,:);
    [pwithin1,pbetween1] = ComputeDiscriminateRBF(D1,D2,bw);
    [pwithin2,pbetween2] = ComputeDiscriminateRBF(D2,D1,bw);

    %compute difference
    Difference1 = ComputeRBFDifference(pwithin1,pbetween1);
    Difference2 = ComputeRBFDifference(pwithin2,pbetween2);
    
    storeit = [storeit; SS, mean([max(Difference1),max(Difference2)] ) ];

end

f = fit(storeit(:,1),storeit(:,2),'exp1');

x = 1:0.1:8;
ab = coeffvalues(f)
z = ab(1).*exp(ab(2).*x);

figure
plot(storeit(:,1),storeit(:,2),'b')
hold on
plot(x,z,'r')
hold off
xlabel('dimensions')
ylabel('max difference')


