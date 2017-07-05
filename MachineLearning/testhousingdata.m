%test machine learning
m = 2.5;
b = 20;
th_act = [b;m]
N = 100;

X = rand(N,1)*10;
X = [ones(length(X),1), X];

Y = X*th_act + rand(N,1)*3;

figure
scatter(X(:,2),Y);



theta = rand(2,1)*10
alpha = 0.01;
numIters = 1000;
thetastore = theta;

%loop through number of iterations
for ii = 1:numIters
    theta = gradientdescent(theta,alpha,Y,X,'batch');
    Error(ii) = costfunction(Y,X,theta);
    thetastore = [thetastore, theta];
end

figure
plot(Error)

disp('numeric theta: ')
theta
disp('numeric error: ')
error1 = costfunction(Y,X,theta)
disp('least squares regression: ')
theta2 = normaleq(X,Y)
disp('ls error: ')
error2 = costfunction(Y,X,theta2)

