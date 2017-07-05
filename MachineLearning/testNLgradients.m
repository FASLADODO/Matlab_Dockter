
cfunc = @(X,Param) Param(1)*sin(X(:,1)) + X(:,2) * Param(2) + exp(X(:,3)*Param(3));

params(1) = 8;
params(2) = 5;
params(3) = 0.1;

N = 100;

X = linspace(1,30,N)';
X = [X,X,X];

%get sample data with noise
Y = cfunc(X,params) + rand(N,1);

%now try to get back params
theta = rand(3,1)*10; %guess
alpha = 0.0001;
numIters = 1000;
thetastore = theta;

%loop through number of iterations
for ii = 1:numIters
    theta = gradientdescentNL(cfunc,theta,Y,X,alpha,'batch');
    Error(ii) = costfunctionNL(Y,X,cfunc,theta);
    thetastore = [thetastore, theta];
end

figure
plot(Error)
title('total error')

disp('numeric theta: ')
theta
disp('actual theta: ')
params

%now plot guess trajectory vs actual
Ybar = cfunc(X,theta);

figure
scatter(X(:,1),Y,'bo');
hold on
scatter(X(:,1),Ybar,'ro');
hold off
title('actual vs GD guess')
legend('actual','guess')