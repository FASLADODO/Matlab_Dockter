% Trying matlabs SVM

%create some fake data
u1 = [8,5];
u2 = [5,4];
c1 = [0.5,0.5];
c2 = [0.5,0.5];
nn = 100;

X0 = [c1(1)*randn(nn,1) + u1(1), c1(2)*randn(nn,1) + u1(2)];
X1 = [c2(1)*randn(nn,1) + u2(1), c2(2)*randn(nn,1) + u2(2)];

%Get working data
X = [X0;X1];
Y = [ones(nn,1)*-1; ones(nn,1)*1]; %(-1/1)

figure
scatter(X0(:,1),X0(:,2),'r.');
hold on
scatter(X1(:,1),X1(:,2),'b.');
hold off
title('orginal data')

w1 = ones(nn,1)/nn;
w2 = ones(nn,1)/nn;
Weights = [w1*0.51; w2*0.49];

SVMModel = fitcsvm(X,Y,'Weights',Weights);

%get models params
W = SVMModel.Beta
B = SVMModel.Bias

%old way
xp = linspace(min(X(:,1)), max(X(:,1)), 100);
yp = - (W(1)*xp + B)/W(2); %THIS IS THE SHIT

%ezplot way
f = @(x,y) B + [x, y]*W;

[P,score] = predict(SVMModel,X);

figure
scatter(X(:,1),X(:,2),20,score(:,1));
hold on
h = ezplot(f, [min(X(:,1)) max(X(:,1)) min(X(:,2)) max(X(:,2))]);
% plot(xp, yp, '-g'); 
hold off
colorbar
colormap cool
title('SVM classify')


%% test  RodsSVM
%also: http://kernelsvm.tripod.com/ (support vector regression)

%create some fake data
u1 = [2,3];
u2 = [3,5];
c1 = [0.5,0.5];
c2 = [0.5,0.5];
nn = 100;

X0 = [c1(1)*randn(nn,1) + u1(1), c1(2)*randn(nn,1) + u1(2)];
X1 = [c2(1)*randn(nn,1) + u2(1), c2(2)*randn(nn,1) + u2(2)];

%Get working data
X = [X0;X1];
Y = [ones(nn,1)*-1; ones(nn,1)*1]; %(-1/1)

%TRy the lagrangian dual
% THIS is thegeneral method
% alpha = zeros(length(Y),1);
% 
% sumtemp = 0;
% for ii = 1:NN
%     for jj = 1:NN
%         sumtemp = sumtemp + alpha(ii)*alpha(jj)*Y(ii)*Y(jj)*X(ii,:)*X(jj,:)' ;
%     end
% end
% 
% W = sum(alpha) - (1/2)*sumtemp;
% %constrain alpha >= 0
% %constrain sum(alpha.*Y) == 0;
%x_i with non-zero alpha_i are called support vectors (SV) 
%solve with quadratic programming


[W,B,SV] = SVM_train(X,Y)

xp = linspace(min(X(:,1)), max(X(:,1)), 100);
yp = - (W(1)*xp + B)/W(2); %THIS IS THE SHIT

%ezplot way
f = @(x,y) B + [x, y]*W';

figure
scatter(X0(:,1),X0(:,2),'r.');
hold on
scatter(X1(:,1),X1(:,2),'b.');
hold on;
h = ezplot(f, [min(X(:,1)) max(X(:,1)) min(X(:,2)) max(X(:,2))]);
% plot(xp, yp, '-g'); 
hold off
title('classify')

[Accuracy,Classify] = SVM_Online(X,Y,W,B);
Accuracy

%% Trying 1d Data

%create some fake data
u1 = [3];
u2 = [5];
c1 = [0.5];
c2 = [0.5];
nn = 100;

X0 = [c1(1)*randn(nn,1) + u1(1)];
X1 = [c2(1)*randn(nn,1) + u2(1)];

%Get working data
X = [X0;X1];
Y = [ones(nn,1)*-1; ones(nn,1)*1]; %(-1/1)


[W,B,SV] = SVM_train(X,Y)

boundary = -B/W %This is the 1D boundary

fig=figure; 
hax=axes;
scatter(X0,ones(nn,1),'r*')
hold on
scatter(X1,ones(nn,1),'b+')
hold on
SP=-B/W; %your point goes here 
line([SP SP],get(hax,'YLim'),'Color',[0 1 0])
hold off

 
[Accuracy,Classify] = SVM_Online(X,Y,W,B);
Accuracy





