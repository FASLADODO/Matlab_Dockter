%Using matlabs built in adaboost sweet

% http://www.mathworks.com/help/stats/compactclassificationtree.view.html
% http://www.mathworks.com/help/stats/ensemble-methods.html
% http://www.mathworks.com/help/stats/fitensemble.html


%fisheriris Doesnt work with AdaBbostM1,Use AdaBoostM2

load fisheriris

%create a single tree
Mdl = fitctree(meas,species);
view(Mdl,'Mode','graph')
pym = predict(Mdl,meas);

%now get the actual data
classes = unique(species);
X = meas;
Y = species;

classes = unique(Y);
Ynum = MapValues(Y,classes,[1,2,3]);

%plot the classes
figure
gscatter3(X(:,1),X(:,2),X(:,3),Ynum);

%now fit the adaboost ensemble tree
ClassTreeEns = fitensemble(X,Y,'AdaBoostM2',100,'Tree');

%compute resub loss
rsLoss = resubLoss(ClassTreeEns,'Mode','Cumulative');

figure
plot(rsLoss);
xlabel('Number of Learning Cycles');
ylabel('Resubstitution Loss');

%try predicting
predMSpec = predict(ClassTreeEns,X); 

%get k fold validation accuracy
cvens = crossval(ClassTreeEns);
L = kfoldLoss(cvens)
 
%% now just some rand data
nn = 100;
ss = 3;

X = [];
Y = [];

%create some random data
s1 = 1;
mu1 = [3,4,2];
s2 = 1;
mu2 = [1,2,1];
x1 = randn(nn,ss).*s1 + repmat(mu1,nn,1);
x2 = randn(nn,ss).*s2 + repmat(mu2,nn,1);

X = [x1;x2];
Y = [repmat({'G1'},nn,1); repmat({'G2'},nn,1) ];

classes = unique(Y);
Ynum = MapValues(Y,classes,[1,2]);

%plot the original classes
figure
gscatter3(X(:,1),X(:,2),X(:,3),Ynum);

%fit adaboost tree
ClassTreeEns = fitensemble(X,Ynum,'AdaBoostM1',20,'Tree');

%get rsloss
rsLoss = resubLoss(ClassTreeEns,'Mode','Cumulative');

figure
plot(rsLoss);
xlabel('Number of Learning Cycles');
ylabel('Resubstitution Loss');


%try predicting 
predY = predict(ClassTreeEns,X);
% corr = strcmp(Y,predY);
corr = Ynum == predY;
sum(corr)/length(corr)

%plot new classes
classes = unique(ClassTreeEns.ClassNames);
predYnum = MapValues(predY,classes,[1,2]);

figure
gscatter3(X(:,1),X(:,2),X(:,3),predYnum);

figure
view(ClassTreeEns.Trained{1},'Mode','graph')


%get rmse
cvens = crossval(ClassTreeEns);
rse = kfoldLoss(cvens)
rsacc = 1- rse

%%

load ionosphere;

ClassTreeEns = fitensemble(X,Y,'AdaBoostM1',100,'Tree');

rsLoss = resubLoss(ClassTreeEns,'Mode','Cumulative');

plot(rsLoss);
xlabel('Number of Learning Cycles');
ylabel('Resubstitution Loss');


%% Make some linear data using 'LSBoost'

nn = 200;
noiz = 0.7;

rangen = [0 25];

params1 = [2.5 10]; %m, b
params2 = [2 10]; %m, b


xm = linspace(rangen(1),rangen(2),nn)';

xm1 = xm + randn(nn,1)*noiz;
xm2 = xm + randn(nn,1)*noiz;

yx1 = xm1.*params1(1) + params1(2) + randn(nn,1)*noiz;
yx2 = xm2.*params2(1) + params2(2) + randn(nn,1)*noiz;


X1 = [xm1, yx1];
X2 = [xm2, yx2];

Data = [X1;X2];
Y = [ones(nn,1)*1;ones(nn,1)*2]; %;ones(nn,1)*3];
classes = unique(Y);

figure
gscatter(Data(:,1),Data(:,2),Y);
hold off
title('Sample Data','FontSize', 12)
xlabel('state1','FontSize', 12)
ylabel('state2','FontSize', 12)
hl = legend('class1','class2','Location','southeast');
set(hl,'FontSize',12);

ClassTreeEns = fitensemble(Data,Y,'LSBoost',100,'Tree');

rsLoss = resubLoss(ClassTreeEns,'Mode','Cumulative');

figure
plot(rsLoss);
xlabel('Number of Learning Cycles');
ylabel('Resubstitution Loss');

predY = predict(ClassTreeEns,Data);

corr = Y == predY;

sum(corr)/length(corr)

figure
gscatter(Data(:,1),Data(:,2),predY);

