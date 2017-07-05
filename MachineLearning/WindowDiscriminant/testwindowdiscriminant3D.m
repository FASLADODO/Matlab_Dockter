%With linear functions

clear all


%%
nn = 1000;
noiz = [0.5, 1, 20];

range1 = [1 30];
range2 = [15 50];

params1 = [-30; 25; 20]; %m, b
params2 = [300; 10; 15]; %m, b


xd1 = linspace(range1(1),range1(2),nn)';
xd2 = linspace(range2(1),range2(2),nn)';

xm1 = [ones(nn,1), xd1 + randn(nn,1)*noiz(1), xd2 + randn(nn,1)*noiz(2)];
xm2 = [ones(nn,1), xd1 + randn(nn,1)*noiz(1), xd2 + randn(nn,1)*noiz(2)];

yx1 = xm1*params1 + randn(nn,1)*noiz(3);
yx2 = xm2*params2 + randn(nn,1)*noiz(3);


X1 = [xm1(:,[2,3]), yx1];
X2 = [xm2(:,[2,3]), yx2];

Data = [X1;X2];
Y = [ones(nn,1)*1;ones(nn,1)*2]; %;ones(nn,1)*3];
classes = unique(Y);

figure
scatter3(X1(:,1),X1(:,2),X1(:,3),'r.')
hold on
scatter3(X2(:,1),X2(:,2),X2(:,3),'b.')

hold off
title('Sample Data','FontSize', 12)
xlabel('state1','FontSize', 12)
ylabel('state2','FontSize', 12)
hl = legend('class1','class2','Location','southeast');
set(hl,'FontSize',12);


%% plot the resulting seperations in 3D
tic
[h] = windowDiscriminantPlot3D(Data, Y, 0.5);
toc

%% Train the models for each
pthresh = 0.5;
kmeanlimit = 7;

if(~exist('kmg'))
    kmg = [];
end
[Models, kmg] = windowDiscriminantTrain3D(Data, Y, pthresh, kmg, kmeanlimit, 'ploton');



%% Log Likelihood ratios with online data

disp('classification accuracy with Rods Discriminant:')
LLR_X1 = windowDiscriminantOnline(X1, Models, classes, pthresh);
LLR_X1(end);
class_okay1 = LLR_X1 >= 0;
Accuracy1 = sum(class_okay1)/length(class_okay1)

LLR_X2 = windowDiscriminantOnline(X2, Models, classes, pthresh);
LLR_X2(end);
class_okay2 = LLR_X2 <= 0;
Accuracy2 = sum(class_okay2)/length(class_okay2)

figure
plot(LLR_X1,'r');
figure
plot(LLR_X2,'b');



