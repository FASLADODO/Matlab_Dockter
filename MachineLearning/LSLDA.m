nn = 100;

mu1 = [ 1,2];
mu2 = [ 5,6];
sigma1 = [1, 0.1; 0.1 ,2.1];
sigma2 = [1, 0.1; 0.1 ,1.5];

data1 = mvnrnd(mu1,sigma1,nn);
data2 = mvnrnd(mu2,sigma2,nn);

X = [data1;data2]; %2D data
Y = [ones(nn,1)*0;ones(nn,1) ]; %0's and 1's


figure
gscatter(X(:,1),X(:,2),Y)

%Least squares LDA
W = inv(X'*X)*X'*Y %get params from least squares approach
B = -0.5; %offset is 0.5

Y_est = X*W + B; %now estimate class

%classify based on sign
class = zeros(length(Y),1);
class(Y_est < 0) = 0;
class(Y_est > 0) = 1;

%get accuracy
acc = Y == class;
sum(acc)/length(acc)

%normal vector
T = [0, -1; 1, 0];
Vec = T*W

%Actual LDA
[Wl, Cl] = LDASimple(X,Y)
[P,Class] = LDAonline(X,Wl,Cl);
Class = MapValues(Class,[1,2],[0,1]);
acc = Y == class;
sum(acc)/length(acc)

[Wl, Cl] = LDASimple(X,Y,P)
[P,Class] = LDAonline(X,Wl,Cl);
Class = MapValues(Class,[1,2],[0,1]);
acc = Y == class;
sum(acc)/length(acc)



figure
plot(Y)
hold on
plot(Y_est)
hold off