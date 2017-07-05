%solve 95% bounding box


nn = 1000;
az = 20;
el = 30;


MU1 = [5 5 5];
SIGMA1 = [3 0.1 0.1; 0.1 3 0.1; 0.1 0.1 3];


X1 = mvnrnd(MU1,SIGMA1,nn);


figure
scatter3(X1(:,1),X1(:,2),X1(:,3),'r.')
hold on
limits = PlotBoundingCube(X1,1.5,[0 1 0])
hold off
xlabel('dim1')
ylabel('dim2')
zlabel('dim3')

%%

classification = BoundingBoxCheck(X1,limits);

figure
scatter3(X1(:,1),X1(:,2),X1(:,3),5,classification)
hold on
PlotBoundingCube(X1,1,[0 1 0]);
hold off
colormap cool
xlabel('dim1')
ylabel('dim2')
zlabel('dim3')

amtwithin = sum(classification)/length(classification)
