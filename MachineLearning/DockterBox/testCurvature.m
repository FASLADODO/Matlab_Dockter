% https://en.wikipedia.org/wiki/Curvature

% k = 1/R should be the recirprocal of radius

xc = 2;
yc = 2;
r = 4;
ang = [0:0.05:2*pi]';

xp=r*cos(ang);
yp=r*sin(ang);
plot(xc+xp,yc+yp);
axis([-4,10,-4,10])

%get derivatives and pad
dX = diff(xp); 
ddX = diff(dX);

dY = diff(yp); 
ddY = diff(dY);

%filter a bit
wsize = 3;
dX = MovingAverage(dX,wsize); 
ddX = MovingAverage(ddX,wsize);
dY = MovingAverage(dY,wsize); 
ddY = MovingAverage(ddY,wsize);


figure
plot(dX)
title('dX')

figure
plot(ddX)
title('ddX')

figure
plot(dY)
title('dY')

figure
plot(ddY)
title('ddY')


dx = dX(1:end-1);
ddx = ddX;
dy = dY(1:end-1);
ddy = ddY;

k = (dx.*ddy - dy.*ddx) ./ ((dx.^2 + dy.^2).^(3/2) );
kact = 1/r

figure
plot(k)
hold on
plot(1:length(k),ones(1,length(k))*kact)
hold off


