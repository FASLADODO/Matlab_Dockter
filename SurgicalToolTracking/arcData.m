% Smallest Arc
steps1 = 0.05;
dist1 = 59.9;

x11 = -dist1:steps1:dist1;
s11 = length(x11)
y11(1:s11) = -dist1;

%%

y12 = -dist1:steps1:dist1;
s12 = length(y12)
x12(1:s12) = -dist1;

%%
yc1 = -82.07;
xc1 = -82.07;
r1 = 143.8;
x13 = -dist1:steps1:dist1;
y13 = yc1 + sqrt(r1*r1 -(x13-xc1).*(x13-xc1));

%%

X1 = [x11,x12,x13];
Y1 = [y11,y12,y13];

scatter(X1,Y1)


%% 
% Middle arc

steps2 = 0.05;
dist2 = 80.0;

x21 = -dist2:steps2:dist2;
s21 = length(x21)
y21(1:s21) = -dist2;

%%

y22 = -dist2:steps2:dist2;
s22 = length(y22)
x22(1:s22) = -dist2;
dist2
%%
yc2 = -90;
xc2 = -90;
r2 = 170.3;
x23 = -dist2:steps2:dist2;
y23 = yc2 + sqrt(r2*r2 -(x23-xc2).*(x23-xc2));

%%

X2 = [x21,x22,x23];
Y2 = [y21,y22,y23];

scatter(X2,Y2)


%% 
% large arc

steps3 = 0.01;
dist3 = 100.0;

x31 = -dist3:steps3:dist3;
s31 = length(x31)
y31(1:s31) = -dist3;

%%

y32 = -dist3:steps3:dist3;
s32 = length(y32)
x32(1:s32) = -dist3;

%%
yc3 = -100;
xc3 = -100;
r3 = 200;
x33 = -dist3:steps3:dist3;
y33 = yc3 + sqrt(r3*r3 -(x33-xc3).*(x33-xc3));

%%

X3 = [x31,x32,x33];
Y3 = [y31,y32,y33];

scatter(X3,Y3)

%%
% plot all
X = [X1,X2,X3];
Y = [Y1,Y2,Y3];

Z = 290; %mm

scatter(X,Y)

