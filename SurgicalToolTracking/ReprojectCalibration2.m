%%
% Test accuracy of point reconstruction


dataDir = [pwd '/myCVdataFolder/']
myExt = '.txt';  

filesInfo = dir(dataDir); 

filesInfo = filesInfo(~[filesInfo.isdir])

fileList = 'DepthCalibration2.txt'

data = load([dataDir fileList]);
size(data)

%%

% Dimensions
width = 640;
height = 480;

xcenter = width / 2.0;
ycenter = height / 2.0;

% Get Middle Point For X and Y from Each Channel
pixX = (data(:,1) + data(:,3)) / 2.0;
pixY = (data(:,2) + data(:,4)) / 2.0;

disp = data(:,5);

realX = data(:,6);
realY = data(:,7);
realZ = data(:,8);

%%

% values for equations
% Z:
% 712.8527   -0.0137   -0.0163   -0.0708    7.1985
% X:
% 0.0016    0.6588
% Y:
% 0.0016   -0.8497

%%

% Az = [24939.2089,-1.046456,-0.00799064,-0.02548524,-14.71241555];

Az = [31992,-1.117];

% indataz = [disp,pixX-xcenter,ycenter-pixY];

indataz = [disp];

ZData = depthCalPower(Az,indataz);


%%

indatax = [ZData,pixX-xcenter];

Ax = [0.0016,0.3643];

XData = XYCal(Ax,indatax);


%%

indatay = [ZData,ycenter-pixY];

Ay = [0.0016,0.0014];

YData = XYCal(Ay,indatay);

%%

error = sqrt((XData - realX).^2 + (YData - realY).^2 + (ZData - realZ).^2);

mean(error)
min(error)
max(error)

%%

% Plot

s = 40;
c1 = [0,0,1];
c2 = [1,0,0];

figure('name','Checking reprojection at 141mm');


scatter3(XData(1:25), YData(1:25), ZData(1:25),s,c1,'x');
axis([-60 60 -60 60 120 160])


hold on

scatter3(realX(1:25),realY(1:25),realZ(1:25),s,c2,'o');

hold off

%%

% Plot

s = 40;
c1 = [0,0,1];
c2 = [1,0,0];

figure('name','Checking reprojection at 181mm');


scatter3(XData(26:50), YData(26:50), ZData(26:50),s,c1,'x');
axis([-75 75 -75 75 120 210])

hold on

scatter3(realX(26:50),realY(26:50),realZ(26:50),s,c2,'o');

hold off

%%

% Plot

s = 40;
c1 = [0,0,1];
c2 = [1,0,0];

figure('name','Checking reprojection at 222mm');


scatter3(XData(51:85), YData(51:85), ZData(51:85),s,c1,'x');
axis([-80 80 -80 80 180 260])


hold on

scatter3(realX(51:85),realY(51:85),realZ(51:85),s,c2,'o');

hold off

%%

% Plot

s = 40;
c1 = [0,0,1];
c2 = [1,0,0];

figure('name','Checking reprojection at 259mm');


scatter3(XData(86:120), YData(86:120), ZData(86:120),s,c1,'x');
axis([-100 100 -100 100 200 300])

hold on

scatter3(realX(86:120),realY(86:120),realZ(86:120),s,c2,'o');

hold off

%%

% Plot

s = 40;
c1 = [0,0,1];
c2 = [1,0,0];

figure('name','Checking reprojection at 299mm');


scatter3(XData(121:155), YData(121:155), ZData(121:155),s,c1,'x');
axis([-150 150 -100 100 250 350])

hold on

scatter3(realX(121:155),realY(121:155),realZ(121:155),s,c2,'o');

hold off