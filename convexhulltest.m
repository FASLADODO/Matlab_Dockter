
d = randn(100,3);
x = d(:,1);
y = d(:,2);
z = d(:,3);
blob = d;
K = convhulln(blob);
figure
scatter3(x,y,z,'b*')
hold on
trisurf(K,blob(:,1),blob(:,2),blob(:,3),'facealpha',0.1)
hold off

%%


pointMatrix = randn(500,3);       %# A set of 20 random 3-D points

[DT,hull] = convHull98Percent(pointMatrix,0.9);

x = pointMatrix(:,1);
y = pointMatrix(:,2);
z = pointMatrix(:,3);  

mask = checkInsideHull(DT,pointMatrix);

figure
scatter3(x,y,z,20,mask,'filled')
hold on
h = convHull3DPlot(DT,hull);
hold off
colormap cool
colorbar

%%

[X,Y,Z] = meshgrid(-4:0.1:4);   %# Create a mesh of coordinates for your volume

simplexIndex = pointLocation(DT,[X(:),Y(:),Z(:)]);
mask = ~isnan(simplexIndex);
mask = reshape(mask,size(X));

figure
scatter3(X(:),Y(:),Z(:),2,mask(:) + 2)
colormap cool
colorbar