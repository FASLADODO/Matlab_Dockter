A = [1,0,0]
B = [1,1,0]

AU = UnitVec(A)
BU = UnitVec(B)

org = [0,0,0];

[angle,axiss,axnorm] = AxisAngle3D(AU,BU)


figure
D = [org;AU];
plot3(D(:,1),D(:,2),D(:,3),'g-')
hold on
D = [org;BU];
plot3(D(:,1),D(:,2),D(:,3),'b-')
hold off
str = sprintf('angle = %f',rad2deg(angle))
title(str)
axis equal