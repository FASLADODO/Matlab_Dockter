clear all

N = 500;

params1(1,1) = 8;
params1(2,1) = 5;

params2(1,1) = 9.5;
params2(2,1) = 5;



X1 = [linspace(-4,20,N)',linspace(1,30,N)'];
X2 = [linspace(5,10,N)',linspace(15,30,N)'];

%get sample data with noise
Y1 = X1*params1 + rand(N,1)*10;
Y2 = X2*params2 + rand(N,1)*10;

figure
scatter(X1(:,1),Y1,'ro')
hold on
scatter(X2(:,1),Y2,'bo')
hold off


[NN,SS] = size(X1);

%now try to get back params
theta1 = rand(SS,1)*10 %guess
theta2 = rand(SS,1)*10 %guess
alpha = 0.001;
numIters = 1000;
thetastore1 = theta1;
thetastore2 = theta2;



%loop through number of iterations
for ii = 1:numIters
    theta1 = gradientdescent(theta1,alpha,Y1,X1,'batch');
    theta2 = gradientdescent(theta2,alpha,Y2,X2,'batch');
    Error1(ii) = costfunction(Y1,X1,theta1);
    Error2(ii) = costfunction(Y2,X2,theta2);
    thetastore1 = [thetastore1, theta1];
    thetastore2 = [thetastore2, theta2];
end

figure
plot(Error1)
title('total error1')

figure
plot(Error2)
title('total error2')

disp('numeric theta: ')
theta1
disp('actual theta: ')
params1

disp('numeric theta: ')
theta2
disp('actual theta: ')
params2

%now plot guess trajectory vs actual
Ybar1 = X1*theta1;

figure
scatter(X1(:,1),Y1,'bo');
hold on
scatter(X1(:,1),Ybar1,'ro');
hold off
title('actual vs GD guess 1')
legend('actual','guess')

Ybar2 = X2*theta2;

figure
scatter(X2(:,1),Y2,'bo');
hold on
scatter(X2(:,1),Ybar2,'ro');
hold off
title('actual vs GD guess 2')
legend('actual','guess')