rand('seed', 912378);

%%
dataDir = [pwd '/myCVdataFolder/']
myExt = '.txt';  

filesInfo = dir(dataDir); 

filesInfo = filesInfo(~[filesInfo.isdir])

fileList = 'calibdatatotal.txt'

data = load([dataDir fileList]);
size(data)


%%

% Dimensions
width = 720;
height = 480;

xcenter = width / 2.0;
ycenter = height / 2.0;

% Get Middle Point For X and Y from Each Channel
X = (data(:,1) + data(:,3)) / 2.0;
Y = (data(:,2) + data(:,4)) / 2.0;

disp = data(:,5);

indata = [disp,X-xcenter,ycenter-Y];

ZData = data(:,8);

%%

% depth vs disparity plot

scatter(data(:,5),data(:,8),15,[0,0,1],'o');

grid on

xlabel('Disparity (pixels)')
ylabel('Z (mm)')

%%

xydata = [X-xcenter,ycenter-Y];

cftool(data(:,8).*xydata(:,1),data(:,6))
cftool(data(:,8).*xydata(:,2),data(:,7))

%%

% find center disparity at each depth

basedepths = [155,175,191,225,250];

dispcenter = [0,0,0,0,0];

keepdex = 0;
found = 0;
for i = 1:length(data(:,1))
   if data(i,6) == -16.15 && data(i,7) == 16.15
       keepdex = i;
       found = 1;
   elseif data(i,6) == 0 && data(i,7) == 16.15
       keepdex = i;
       found = 1;
   elseif data(i,6) == 16.15 && data(i,7) == 16.15
       keepdex = i;
       found = 1;
   elseif data(i,6) == -16.15 && data(i,7) == 0
       keepdex = i;
       found = 1;
   elseif data(i,6) == 16.15 && data(i,7) == 0
       keepdex = i;
       found = 1;
   elseif data(i,6) == -16.15 && data(i,7) == -16.15
       keepdex = i;
       found = 1;
   elseif data(i,6) == 0 && data(i,7) == -16.15
       keepdex = i;
       found = 1;
   elseif data(i,6) == 16.15 && data(i,7) == -16.15
       keepdex = i;
       found = 1;
   end
   
   if found == 1
      if(data(keepdex,8) == basedepths(1))
          dispcenter(1) = dispcenter(1)+data(keepdex,5);
      elseif(data(keepdex,8) == basedepths(2))
          dispcenter(2) = dispcenter(2)+data(keepdex,5);
      elseif(data(keepdex,8) == basedepths(3))
          dispcenter(3) = dispcenter(3)+data(keepdex,5);
      elseif(data(keepdex,8) == basedepths(4))
          dispcenter(4) = dispcenter(4)+data(keepdex,5);
      elseif(data(keepdex,8) == basedepths(5))
          dispcenter(5) = dispcenter(5)+data(keepdex,5);
      end
   end
   
   keepdex = 0;
   found = 0;
    
end

totdisp = dispcenter/8

dispfit = zeros(length(data(:,1)),1);

for j = 1:length(data(:,1))
   if data(j,8) == basedepths(1)
       dispfit(j,1) = totdisp(1);
   elseif data(j,8) == basedepths(2)
       dispfit(j,1) = totdisp(2);    
   elseif data(j,8) == basedepths(3)
       dispfit(j,1) = totdisp(3);
   elseif data(j,8) == basedepths(4)
       dispfit(j,1) = totdisp(4);  
   elseif data(j,8) == basedepths(5)
       dispfit(j,1) = totdisp(5);     
   end
    
end

%%
[r1,c,v] = find(ZData == 155);
[r2,c,v] = find(ZData == 175);
[r3,c,v] = find(ZData == 191);
[r4,c,v] = find(ZData == 225);
[r5,c,v] = find(ZData == 250);

scatter3(X(r1,1),Y(r1,1),data(r1,5)-totdisp(1))

xlabel('x [pixels]')
ylabel('y [pixels]')
zlabel('disparity [pixel delta]')

%%

% random fits

N = 99928;
%rseed(N);

fcoeff = 16;
M = 100000; % #loops

rc = (rand(fcoeff,M)-0.5).*1.*2;

min(rc(1,:))
max(rc(1,:))

xydata = [X-xcenter,ycenter-Y];

fitdata = data(:,5)-dispfit;

min(fitdata)
max(fitdata)
length(fitdata)

min(xydata)
max(xydata)
length(xydata)


%%
AP = [-0.2782,-0.01216,-0.002626,3.312e-05,6.151e-05,9.884e-06,5.155e-07,9.715e-08,5.644e-07,6.742e-09];
   
[params,R,J,covB,MSE] = nlinfit(xydata,fitdata,@xyoffsetfit,AP);

params

%%

% r1
% Linear model Poly33:
%      f(x,y) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 + p21*x^2*y 
%                     + p12*x*y^2 + p03*y^3
% Coefficients (with 95% confidence bounds):
%        p00 =     -0.5779  (-1.13, -0.02628)
%        p10 =    -0.01343  (-0.01797, -0.008886)
%        p01 =  -0.0005261  (-0.009764, 0.008712)
%        p20 =   4.188e-05  (3.088e-05, 5.288e-05)
%        p11 =   6.979e-05  (5.234e-05, 8.723e-05)
%        p02 =   2.211e-05  (-1.762e-05, 6.185e-05)
%        p30 =   5.652e-07  (4.901e-07, 6.403e-07)
%        p21 =    7.89e-08  (-3.548e-08, 1.933e-07)
%        p12 =   5.435e-07  (3.248e-07, 7.622e-07)
%        p03 =  -4.225e-08  (-5.894e-07, 5.049e-07)
% 
% Goodness of fit:
%   SSE: 26.49
%   R-square: 0.9701
%   Adjusted R-square: 0.9621
%   RMSE: 0.8827

%%

% r2 Linear model Poly33:
%      f(x,y) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 + p21*x^2*y 
%                     + p12*x*y^2 + p03*y^3
% Coefficients (with 95% confidence bounds):
%        p00 =      -0.224  (-0.5344, 0.08639)
%        p10 =   -0.009104  (-0.01195, -0.00626)
%        p01 =   -0.005672  (-0.01143, 8.512e-05)
%        p20 =   3.245e-05  (2.478e-05, 4.013e-05)
%        p11 =   5.168e-05  (3.95e-05, 6.386e-05)
%        p02 =   3.394e-05  (6.03e-06, 6.186e-05)
%        p30 =   4.016e-07  (3.434e-07, 4.598e-07)
%        p21 =   2.075e-07  (1.19e-07, 2.961e-07)
%        p12 =   5.967e-07  (4.28e-07, 7.655e-07)
%        p03 =   5.441e-08  (-3.674e-07, 4.762e-07)
% 
% Goodness of fit:
%   SSE: 8.349
%   R-square: 0.9681
%   Adjusted R-square: 0.9597
%   RMSE: 0.4955
%%
% r3

% Linear model Poly33:
%      f(x,y) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 + p21*x^2*y 
%                     + p12*x*y^2 + p03*y^3
% Coefficients (with 95% confidence bounds):
%        p00 =     -0.4409  (-0.7937, -0.08808)
%        p10 =    -0.01183  (-0.01544, -0.008221)
%        p01 =   -0.003259  (-0.008044, 0.001527)
%        p20 =   3.319e-05  (2.258e-05, 4.38e-05)
%        p11 =   6.839e-05  (5.649e-05, 8.029e-05)
%        p02 =   1.745e-05  (-7.682e-07, 3.566e-05)
%        p30 =   5.121e-07  (4.241e-07, 6.001e-07)
%        p21 =   1.315e-07  (3.645e-08, 2.266e-07)
%        p12 =   5.539e-07  (4.296e-07, 6.782e-07)
%        p03 =  -2.893e-08  (-2.274e-07, 1.695e-07)
% 
% Goodness of fit:
%   SSE: 25.16
%   R-square: 0.9388
%   Adjusted R-square: 0.9282
%   RMSE: 0.6956

%%
%r4

% Linear model Poly33:
%      f(x,y) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 + p21*x^2*y 
%                     + p12*x*y^2 + p03*y^3
% Coefficients (with 95% confidence bounds):
%        p00 =     -0.1074  (-0.4737, 0.2589)
%        p10 =   -0.008655  (-0.0122, -0.005114)
%        p01 =   0.0003411  (-0.005448, 0.00613)
%        p20 =   3.264e-05  (2.276e-05, 4.253e-05)
%        p11 =   5.449e-05  (4.071e-05, 6.827e-05)
%        p02 =  -1.504e-05  (-4.102e-05, 1.093e-05)
%        p30 =   4.497e-07  (3.726e-07, 5.268e-07)
%        p21 =   1.043e-08  (-9.299e-08, 1.138e-07)
%        p12 =   4.364e-07  (2.701e-07, 6.026e-07)
%        p03 =  -1.023e-07  (-4.265e-07, 2.219e-07)
% 
% Goodness of fit:
%   SSE: 43.07
%   R-square: 0.9129
%   Adjusted R-square: 0.901
%   RMSE: 0.8078

%%

% r5
% Linear model Poly33:
%      f(x,y) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 + p21*x^2*y 
%                     + p12*x*y^2 + p03*y^3
% Coefficients (with 95% confidence bounds):
%        p00 =     -0.1366  (-0.4764, 0.2033)
%        p10 =    -0.01478  (-0.0185, -0.01106)
%        p01 =    -0.00335  (-0.007935, 0.001235)
%        p20 =   2.168e-05  (1.002e-05, 3.335e-05)
%        p11 =   5.962e-05  (4.705e-05, 7.219e-05)
%        p02 =   1.214e-05  (-5.865e-06, 3.014e-05)
%        p30 =   5.834e-07  (4.814e-07, 6.854e-07)
%        p21 =   1.061e-07  (2.1e-10, 2.121e-07)
%        p12 =    6.25e-07  (4.946e-07, 7.553e-07)
%        p03 =   5.992e-08  (-1.306e-07, 2.504e-07)
% 
% Goodness of fit:
%   SSE: 65.92
%   R-square: 0.8644
%   Adjusted R-square: 0.8505
%   RMSE: 0.8655

%%

% total

% Linear model Poly33:
%      f(x,y) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 + p21*x^2*y 
%                     + p12*x*y^2 + p03*y^3
% Coefficients (with 95% confidence bounds):
%        p00 =     -0.2782  (-0.4498, -0.1066)
%        p10 =    -0.01216  (-0.01378, -0.01053)
%        p01 =   -0.002626  (-0.005051, -0.0002016)
%        p20 =   3.312e-05  (2.844e-05, 3.78e-05)
%        p11 =   6.151e-05  (5.525e-05, 6.777e-05)
%        p02 =   9.884e-06  (-4.667e-07, 2.023e-05)
%        p30 =   5.155e-07  (4.814e-07, 5.495e-07)
%        p21 =   9.715e-08  (5.016e-08, 1.441e-07)
%        p12 =   5.644e-07  (4.948e-07, 6.34e-07)
%        p03 =   6.742e-09  (-1.042e-07, 1.177e-07)
% 
% Goodness of fit:
%   SSE: 202.9
%   R-square: 0.9206
%   Adjusted R-square: 0.9183
%   RMSE: 0.8039
%%

% average

p00 =     -0.5779 - 0.224 - 0.4409 -0.1074 - 0.1366;
p10 =    -0.01343 - 0.009104 - 0.01183-0.008655 - 0.01478;
p01 =  -0.0005261 - 0.005672 - 0.003259 + 0.0003411 - 0.00335;
p20 =   4.188e-05 + 3.245e-05 + 3.319e-05 + 3.264e-05 + 2.168e-05;
p11 =   6.979e-05 + 5.168e-05 + 6.839e-05 + 5.449e-05 + 5.962e-05;
p02 =   2.211e-05 + 3.394e-05 + 1.745e-05 -1.504e-05 + 1.214e-05;
p30 =   5.652e-07 + 4.016e-07 + 5.121e-07 + 4.497e-07 + 5.834e-07;
p21 =    7.89e-08 + 2.075e-07 + 1.315e-07 + 1.043e-08 + 1.061e-07;
p12 =   5.435e-07 + 5.967e-07 + 5.539e-07 + 4.364e-07 + 6.25e-07;
p03 =  -4.225e-08 + 5.441e-08 - 2.893e-08 + -1.023e-07 + 5.992e-08; 

p00 = p00 / 5
p10 = p10 / 5
p01 = p01 / 5
p20 = p20 / 5
p11 = p11 / 5
p02 = p02 / 5
p30 = p30 / 5
p21 = p21 / 5
p12 = p12 / 5
p03 = p03 / 5



%%

testdata = [X-xcenter,ycenter-Y];

fitdata = data(:,5)-dispfit;

cftool(testdata(r1,1),testdata(r1,2),fitdata(r1,1))
%%
cftool(testdata(r2,1),testdata(r2,2),fitdata(r2,1))
%%
cftool(testdata(r3,1),testdata(r3,2),fitdata(r3,1))

%%
cftool(testdata(r4,1),testdata(r4,2),fitdata(r4,1))

%%
cftool(testdata(r5,1),testdata(r5,2),fitdata(r5,1))

%%
cftool(testdata(:,1),testdata(:,2),fitdata(:,1))


%%

% using cftool

s = 30;

c1 = [0,0,1];
c2 = [1,0,0];
c3 = [0,1,0];
c4 = [1,1,0];
c5 = [1,0,1];

[r1,c,v] = find(ZData == 155);
[r2,c,v] = find(ZData == 175);
[r3,c,v] = find(ZData == 191);
[r4,c,v] = find(ZData == 225);
[r5,c,v] = find(ZData == 250);

indivdepths = [155,175,191,225,250];

testdatar = [X-xcenter,ycenter-Y];

APx = [0.00155,-2.495];

APy = [0.001547,2.674];

AP2 = [-0.2974,-0.0116,-0.0025,3.2368e-05,6.0794e-05,1.4120e-05,5.0240e-07,1.0689e-07,5.5110e-07,-1.1830e-08];

AP = [-0.2782,-0.01216,-0.002626,3.312e-05,6.151e-05,9.884e-06,5.155e-07,9.715e-08,5.644e-07,6.742e-09];

APz = [7.615,0.05419];

offsetter = xyoffsetfit(AP,[0,0])

fitdata = 54-offsetter

calcdepth = depthCalDaVinci(APz,fitdata)

calcx = XYCalDaVinci(APx,[calcdepth,0])
calcy = XYCalDaVinci(APy,[calcdepth,0])


%%

offsetter = xyoffsetfit(AP,testdatar);

fitdata = data(:,5)-offsetter;

testthisshit = fitdata - dispfit;

% cftool(testdatar(:,1),testdatar(:,2),testthisshit);

calcdepth = depthCalDaVinci(APz,fitdata);

calcx = XYCalDaVinci(APx,[calcdepth,testdatar(:,1)]);
calcy = XYCalDaVinci(APy,[calcdepth,testdatar(:,2)]);

figure('name','reprojected vs model 150');

scatter3(calcx(r1,1),calcy(r1,1),calcdepth(r1,1),s,c1,'x');

hold on

scatter3(data(r1,6),data(r1,7),data(r1,8),s,c2,'o');

hold off
xlabel('x [world]')
ylabel('y [world]')
zlabel('Z [world]')
 legend('data','known')

%%

figure('name','reprojected vs model 175');

scatter3(calcx(r2,1),calcy(r2,1),calcdepth(r2,1),s,c1,'x');

hold on

scatter3(data(r2,6),data(r2,7),data(r2,8),s,c2,'o');

hold off
xlabel('x [world]')
ylabel('y [world]')
zlabel('Z [world]')
legend('data','known')

%%

figure('name','reprojected vs model 191');

scatter3(calcx(r3,1),calcy(r3,1),calcdepth(r3,1),s,c1,'x');

hold on

scatter3(data(r3,6),data(r3,7),data(r3,8),s,c2,'o');

hold off
xlabel('x [world]')
ylabel('y [world]')
zlabel('Z [world]')
legend('data','known')

%%

figure('name','reprojected vs model 225');

scatter3(calcx(r4,1),calcy(r4,1),calcdepth(r4,1),s,c1,'x');

hold on

scatter3(data(r4,6),data(r4,7),data(r4,8),s,c2,'o');

hold off
xlabel('x [world]')
ylabel('y [world]')
zlabel('Z [world]')
legend('data','known')

%%

figure('name','reprojected vs model 250');

scatter3(calcx(r5,1),calcy(r5,1),calcdepth(r5,1),s,c1,'x');

hold on

scatter3(data(r5,6),data(r5,7),data(r5,8),s,c2,'o');

hold off
xlabel('x [world]')
ylabel('y [world]')
zlabel('Z [world]')
legend('data','known')

%%

errorxy = sqrt((calcx(:,1)-data(:,6)).^2 + (calcy(:,1)-data(:,7)).^2);

error = sqrt((calcx(:,1)-data(:,6)).^2 + (calcy(:,1)-data(:,7)).^2 + (calcdepth(:,1)-data(:,8)).^2);

mean(error)

std(error)
std(errorxy)

prctile(error,95)

%%

indivdepths = [155,175,191,225,250];

s = 10;

c1 = [0,0,1];
c2 = [1,0,0];
c3 = [0,1,0];
c4 = [1,1,0];
c5 = [1,0,1];

figure('name','X vs Y vs disp');

scatter3(testdatar(r1,1),testdatar(r1,2),fitdata(r1,1),s,c1);

hold on
%
scatter3(testdatar(r2,1),testdatar(r2,2),fitdata(r2,1),s,c2);

hold on
scatter3(testdatar(r3,1),testdatar(r3,2),fitdata(r3,1),s,c3);

hold on
scatter3(testdatar(r4,1),testdatar(r4,2),fitdata(r4,1),s,c4);

hold on
scatter3(testdatar(r5,1),testdatar(r5,2),fitdata(r5,1),s,c5);

hold off
xlabel('x [pixels]')
ylabel('y [pixels]')
zlabel('disparity [pixel delta]')
Txt = [];
for i =1:length(indivdepths)
    Txt{i} = [num2str(indivdepths(i)) ' mm'];
end
    
legend(Txt);

%%

%first depth
[r1,c,v] = find(ZData == 155);

figure('name','X vs disp vs Depth 155');

scatter3(data(r1,5),X(r1),data(r1,8))

%% second depth
[r2,c,v] = find(ZData == 175);

figure('name','X vs disp vs Depth 175');

scatter3(data(r2,5),X(r2),data(r2,8))

%% third depth
[r3,c,v] = find(ZData == 191);

figure('name','X vs disp vs Depth 191');

scatter3(data(r3,5),X(r3),data(r3,8))

%% fourth depth

[r4,c,v] = find(ZData == 225);

figure('name','X vs disp vs Depth 225');

scatter3(data(r4,5),X(r4),data(r4,8))

%% fifth depth

[r5,c,v] = find(ZData == 250);

figure('name','X vs disp vs Depth 250');

scatter3(data(r5,5),X(r5),data(r5,8))
