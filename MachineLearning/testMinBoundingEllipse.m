
%create some rotated gaussian data
cv = [0.8,0.3;
    0.3,0.8];
mu = [3,4];
points = mvnrnd(mu,cv,10000);

%get the minimum bounding ellipse
[A, centroid] = MinBoundEllipse_Compute(points',0.1);

A
centroid


%plot the data and centroid
figure
plot(points(:,1),points(:,2),'.');
hold on;
h = MinBoundEllipse_Plot(A,centroid,50);
axis square;

%Test how many points are inside the ellipse
P = MinBoundEllipse_Test(A,centroid,points)