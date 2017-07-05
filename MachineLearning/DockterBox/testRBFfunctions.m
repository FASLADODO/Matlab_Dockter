
nn = 50;
t = linspace(1,2,nn)';

X = t.^7;



%get derivatives and pad
dX = diff(X); 
ddX = diff(dX);
dddX = diff(ddX);


%add some noise
dX = dX + randn(nn-1,1)*0.1;
ddX = ddX + randn(nn-2,1)*0.1;

%filter a bit
wsize = 5;
dX = MovingAverage(dX,wsize); 
ddX = MovingAverage(ddX,wsize);
dddX = MovingAverage(dddX,wsize);

figure
scatter(t,X)
title('X')

figure
scatter(t(1:end-1),dX)
title('dX')

figure
scatter(t(1:end-2),ddX)
title('ddX')

figure
scatter(t(1:end-3),dddX)
title('dddX')

data = [dX(1:end-1,:),ddX];
gamma =2;
Z = RadialBasisFunction(data,gamma);

figure
scatter3(data(:,1),data(:,2),Z)

figure
scatter(t(1:end-2),Z)
title('RBF t^5')
ylim([0,1])

mean(Z)

%% Test ordering

nn = 10;
x1 = linspace(1,10,nn)';
x2 = linspace(1,20,nn)';

x2 = x2 + randn(nn,1);

[p,S] = polyfit(x1,x2,2)
yf = polyval(p,x1);
[r2 rmse] = rsquare(x2,yf)

XR = [ x1,x2];
ZR = RadialBasisFunction(XR,gamma);

figure
scatter(XR(:,1),XR(:,2))
hold on
plot(x1,yf)
hold off
title('XR')

figure
scatter3(XR(:,1),XR(:,2),ZR)



