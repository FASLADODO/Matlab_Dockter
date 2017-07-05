%create some fake data
u1 = [2,3];
u2 = [2,5];
c1 = [0.1,0.9];
c2 = [0.9,0.1];
nn = 300;

X0 = [c1(1)*randn(nn,1) + u1(1), c1(2)*randn(nn,1) + u1(2)];
X1 = [c2(1)*randn(nn,1) + u2(1), c2(2)*randn(nn,1) + u2(2)];


%Get working data
X = [X0,X1];
X = sort(X,1);
noise = randn(nn,1)*0.5;

%Just random numbers from a model
Param_True = [1;2;-7;1]; 

%get Y
Y = X*Param_True + noise;

%plot Y
figure
scatter(1:nn,Y)
xlabel('samples')
ylabel('Y(x)')
title('sample data')

%back convert params for science
param_test = inv(X'*X)*X'*Y

%Do bayesian information criteria
[FullErrors, AvgErrors, paramvars] = BayesianInformationCriteria(X,Y,{'x1','x2','x3','x4'},'sample data','TLS');

