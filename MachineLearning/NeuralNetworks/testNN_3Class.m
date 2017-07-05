%3D machine learning
clear all



%% 3D

nn = 50;
az = 20;
el = 30;


MU1 = [0 0 0];
SIGMA1 = [1 0.1 0.1; 0.1 1 0.1; 0.1 0.1 1];

MU2 = [0 5 6];
SIGMA2 = [1 0.1 0.1; 0.1 1 0.1; 0.1 0.1 1];

MU3 = [7 4 0];
SIGMA3 = [1 0.1 0.1; 0.1 1 0.1; 0.1 0.1 1];

X1 = mvnrnd(MU1,SIGMA1,nn);
X2 = mvnrnd(MU2,SIGMA2,nn);
X3 = mvnrnd(MU3,SIGMA3,nn);

X = [X1;X2;X3];
Y = [ones(nn,1)*1; ones(nn,1)*2; ones(nn,1)*3];
classes = unique(Y);

figure
gscatter3(X(:,1),X(:,2),X(:,3),Y)

title('Sample Data','FontSize', 12)
xlabel('state1','FontSize', 12)
ylabel('state2','FontSize', 12)



%% Neural Nets

layers = 1;
nodes = 3;

%we have to format because too many classes
[YBinary,cslist,mapping] = NNFormatOutput(Y);

[NN] = NNInitialize(X,YBinary,layers,nodes);

alpha = 0.5; %learning rate (0 - 1)
max_epoch = 1000;
accuracy_target = 0.01;

[NN,MSE] = NNTrain(X,YBinary,NN,alpha,max_epoch,accuracy_target);

%% Check classifiying

%online NN classify
Yestb = NNOnline(NN, X);

%convert back to our class labels
Yest = NNUnformatOutput(Yestb,cslist,mapping);
    
figure
gscatter3(X(:,1),X(:,2),X(:,3),Yest)













