%create data matrix
nn = 1000;
d = 3;
mu1 = [0 1 0];
scale1 = [4 4 4];
mu2 = [0 1 0];
scale2 = [2 2 2];
mu3 = [0 1 0];
scale3 = [1 1 1];

%create data and plot
dat1 = [randn(nn,1)*(scale1(1)/2) + mu1(1), randn(nn,1)*(scale1(2)/2) + mu1(2), randn(nn,1)*(scale1(3)/2) + mu1(3)];
dat2 = [randn(nn,1)*(scale2(1)/2) + mu2(1), randn(nn,1)*(scale2(2)/2) + mu2(2), randn(nn,1)*(scale2(3)/2) + mu2(3)];
dat3 = [randn(nn,1)*(scale3(1)/2) + mu3(1), randn(nn,1)*(scale3(2)/2) + mu3(2), randn(nn,1)*(scale3(3)/2) + mu3(3)];
dat = [dat1;dat2;dat3];

spherecolor = [1 1 0];


%class 1
figure
[radius, center] = bounding_sphere(dat1,0.95)
[X,Y,Z] = sphereObject(radius,center);
scatter3(dat1(:,1),dat1(:,2),dat1(:,3),'b.')
hold on
m=mesh(X,Y,Z);
set(m,'facecolor','none')
colormap(spherecolor)
hold off
title('data 1')
xlabel('pos')
ylabel('velocity')
zlabel('acceleration')
legend('data','bounding')
axis([-7 7 -7 7 -7 7])

%check if actually 95%
poo = [];
daton = dat1; %get data
datc = daton - center(ones(size(daton,1),1),:); %centerthatshit
mag = sqrt(sum(datc.^2, 2)); %distance
radz = ones(length(datc),1)*radius; %vector of radius for comparing
poo = mag <= radz; %compare it
acc = sum(poo)/length(poo) %accuracy


%class 2
figure
[radius, center] = bounding_sphere(dat2,0.95)
[X,Y,Z] = sphereObject(radius,center);
scatter3(dat2(:,1),dat2(:,2),dat2(:,3),'r.')
hold on
m=mesh(X,Y,Z);
set(m,'facecolor','none')
colormap(spherecolor)
hold off
title('data 2')
xlabel('pos')
ylabel('velocity')
zlabel('acceleration')
legend('data','bounding')
axis([-7 7 -7 7 -7 7])

%check if actually 95%
poo = [];
daton = dat2; %get data
datc = daton - center(ones(size(daton,1),1),:); %centerthatshit
mag = sqrt(sum(datc.^2, 2)); %distance
radz = ones(length(datc),1)*radius; %vector of radius for comparing
poo = mag <= radz; %compare it
acc = sum(poo)/length(poo) %accuracy


%class 3
figure
[radius, center] = bounding_sphere(dat3,0.95)
[X,Y,Z] = sphereObject(radius,center);
scatter3(dat3(:,1),dat3(:,2),dat3(:,3),'k.')
hold on
m=mesh(X,Y,Z);
set(m,'facecolor','none')
colormap(spherecolor)
hold off
title('data 3')
xlabel('pos')
ylabel('velocity')
zlabel('acceleration')
legend('data','bounding')
axis([-7 7 -7 7 -7 7])

%check if actually 95%
poo = [];
daton = dat3; %get data
datc = daton - center(ones(size(daton,1),1),:); %centerthatshit
mag = sqrt(sum(datc.^2, 2)); %distance
radz = ones(length(datc),1)*radius; %vector of radius for comparing
poo = mag <= radz; %compare it
acc = sum(poo)/length(poo) %accuracy


