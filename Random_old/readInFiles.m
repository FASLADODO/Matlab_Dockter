%% Read in a bunch of files from a given directory

% target dir where data lives (assume only data files have extension myExt)
dataDir = [pwd '/myCVdataFolder/']
myExt = '.txt';  % assume data files have this extension

% array of fileInof structs
filesInfo = dir(dataDir); 

%remove directories
filesInfo = filesInfo(~[filesInfo.isdir])

% get list of filenames to open

fileList = {}; 

% Read in bunches of files
% for k=1:length(filesInfo)
%     % if it matches extension, save it.
%     if strfind(filesInfo(k).name, myExt)
%         fileList{end+1} = [  filesInfo(k).name]
%     end
% end

%use one file
fileList{1} = 'ToolDataDv200 (198).txt'
fileList{2} = 'ToolDataDv214 (212).txt'
fileList{3} = 'ToolDataDv227 (225).txt'
fileList{4} = 'ToolDataDv245 (246).txt'
fileList{5} = 'ToolDataDv250 (248).txt'
fileList{6} = 'ToolDataDv252 (251).txt'
fileList{7} = 'ToolDataDv280 (278).txt'

%%

data1 = load([dataDir fileList{1}]);
data2 = load([dataDir fileList{2}]);
data3 = load([dataDir fileList{3}]);
data4 = load([dataDir fileList{4}]);
data5 = load([dataDir fileList{5}]);
data6 = load([dataDir fileList{6}]);
data7 = load([dataDir fileList{7}]);

data = data7;
size(data)

%compute overall metrics
fprintf('Overall Means \n');

capdiff = abs(data(:,15) - data(:,14));
trackRtime = abs(data(:,17) - data(:,16));
trackLtime = abs(data(:,19) - data(:,18));
depthTime = abs(data(:,21) - data(:,20));
totalTime = abs(data(:,21) - data(:,14));

% disp('Capture Time Difference Mean')
meancap = mean(capdiff)

% disp('Tracking Right Time Mean')
trackright = mean(trackRtime)

% disp('Tracking Left Time Mean')
trackleft = mean(trackLtime)

% disp('Depth Time Mean')
depthmean = mean(depthTime)

% disp('Total Time Right Mean')
totalright = mean(trackRtime+depthTime)

% disp('Total Time Left Mean')
totalleft = mean(trackLtime+depthTime)

% disp('Total Time Log Mean')
totallog = mean(totalTime)

% disp('Implied FPS')
impfps = 1 / ((max(totalright,totalleft) / 1000))


impfps2 = 1 / ((max(trackright,trackleft) / 1000))


mean(data(:,4))

%%

min(totalTime)
max(totalTime)

%%

datasize = length(data(:,2))


xavg = mean(data(:,2))
yavg = mean(data(:,3))
zavg = mean(data(:,4))

Xdata = data(:,2)';
Ydata = data(:,3)';
Zdata = data(:,4)';

xdist = zeros(1,datasize);
ydist = zeros(1,datasize);
zdist = zeros(1,datasize);
totaldist = zeros(1,datasize);

for j = 1:datasize
    
    totaldist(j) = sqrt((xavg-Xdata(j))^2 + (yavg-Ydata(j))^2 + (zavg-Zdata(j))^2);

    xdist(j) = abs(xavg - Xdata(j));
    
    ydist(j) = abs(yavg - Ydata(j));
    
    zdist(j) = abs(zavg - Zdata(j));
    
end

xdev = mean(xdist)

ydev = mean(ydist)

zdev = mean(zdist)

totaldev = mean(totaldist)

s1 = 10;
s2 = 90;
c1 = [0,0,1];
c2 = [1,0,0];

figure('name','Deviation');

scatter3(Xdata,Ydata,Zdata,s1,c1,'o');

axis([xavg-5,xavg+5,yavg-5,yavg+5,zavg-5,zavg+5])

hold on

scatter3(xavg,yavg,zavg,s2,c2,'x');

hold off

xlabel('x [mm]')
ylabel('y [mm]')
zlabel('z [mm]')
legend('data','point')

%%
[val,ind] = max(data(:,4))
min(data(:,4))

%%

% Smallest Arc
steps1 = 0.05;
dist1 = 59.9;

x11 = -dist1:steps1:dist1;
s11 = length(x11)
y11(1:s11) = -dist1;

y12 = -dist1:steps1:dist1;
s12 = length(y12)
x12(1:s12) = -dist1;

yc1 = -82.07;
xc1 = -82.07;
r1 = 143.8;
x13 = -dist1:steps1:dist1;
y13 = yc1 + sqrt(r1*r1 -(x13-xc1).*(x13-xc1));


X1 = [x11,x12,x13];
Y1 = [y11,y12,y13];


stotal1 = length(X1)

Z11(1:stotal1) = 198.55;
Z12(1:stotal1) = 212.12;
Z13(1:stotal1) = 225.46;
Z14(1:stotal1) = 245.97;
Z15(1:stotal1) = 248.85;
Z16(1:stotal1) = 251.03;
Z17(1:stotal1) = 285.69;

length(Z12)
length(X1)
length(Y1)

%%
% Middle arc

steps2 = 0.05;
dist2 = 80.0;

x21 = -dist2:steps2:dist2;
s21 = length(x21)
y21(1:s21) = -dist2;


y22 = -dist2:steps2:dist2;
s22 = length(y22)
x22(1:s22) = -dist2;
dist2

yc2 = -90;
xc2 = -90;
r2 = 170.3;
x23 = -dist2:steps2:dist2;
y23 = yc2 + sqrt(r2*r2 -(x23-xc2).*(x23-xc2));


X2 = [x21,x22,x23];
Y2 = [y21,y22,y23];

stotal2 = length(X2)

Z21(1:stotal2) = 289.77;
Z22(1:stotal2) = 305.851;
Z23(1:stotal2) = 329.43;
Z24(1:stotal2) = 280.02;
Z25(1:stotal2) = 273.4;

length(Z24)
length(X2)
length(Y2)

%% 
% large arc

steps3 = 0.05;
dist3 = 100.0;

x31 = -dist3:steps3:dist3;
s31 = length(x31)
y31(1:s31) = -dist3;


y32 = -dist3:steps3:dist3;
s32 = length(y32)
x32(1:s32) = -dist3;

yc3 = -100;
xc3 = -100;
r3 = 200;
x33 = -dist3:steps3:dist3;
y33 = yc3 + sqrt(r3*r3 -(x33-xc3).*(x33-xc3));


X3 = [x31,x32,x33];
Y3 = [y31,y32,y33];

stotal3 = length(X3)

Z31(1:stotal3) = 308.636;

length(Z31)
length(X3)
length(Y3)

%%

%Plot arcs

s = 10;
c1 = [0,0,1];
c2 = [1,0,0];

figure('name','All Arcs');

scatter3(X1,Y1,Z1,s,c2);

axis([-120,120,-120,120,200,300])

hold on

scatter3(X2,Y2,Z2,s,c2);

hold on

scatter3(X3,Y3,Z3,s,c2);

hold off

%%

%Plot single arc

Xdata = data(:,2)';
Ydata = data(:,3)';
Zdata = data(:,4)';

s1 = 20;
s2 = 5
c1 = [0,0,1];
c2 = [1,0,0];

figure('name','Small Arc');

scatter3(data(:,2),data(:,3),data(:,4),s1,c1,'x');

axis([-100,100,-100,100,200,300])

hold on

scatter3(X1,Y1,Z11,s2,c2,'o');

hold off
xlabel('x [mm]')
ylabel('y [mm]')
zlabel('Z [mm]')
legend('data','model')


%%

hold on

scatter3(X2,Y2,Z2,s,c2);

hold on

scatter3(X3,Y3,Z3,s,c2);

hold off

%%

%avgerr for small arc 4.4215

%avgerr for all 3 = 21.2 mm

%%
%single arc

% Z11(1:stotal1) = 198.55;
% Z12(1:stotal1) = 212.12;
% Z13(1:stotal1) = 225.46;
% Z14(1:stotal1) = 245.97;
% Z15(1:stotal1) = 248.85;
% Z16(1:stotal1) = 251.03;
% Z17(1:stotal1) = 285.69;

data = data7;

mindist = 1000;
mindistx = 1000;
mindisty = 1000;
mindistz = 1000;

saveindex = 1;
saveindexy = 1;
saveindexx = 1;
saveindexz = 1;

XModel = X1;
YModel = Y1;
ZModel = Z17;

Xdata = data(:,2)';
Ydata = data(:,3)';
Zdata = data(:,4)';

lengthdata = length(Xdata)

lengthmodel = length(XModel)


Xsave = zeros(1,lengthdata);
Ysave = zeros(1,lengthdata);
Zsave = zeros(1,lengthdata);

error = zeros(1,lengthdata);
errorx = zeros(1,lengthdata);
errory = zeros(1,lengthdata);
errorz = zeros(1,lengthdata);

percentcnt = 0;

for i = 1:lengthdata
    for j = 1:lengthmodel
        dist = sqrt((Xdata(i)-XModel(j))^2 + (Ydata(i)-YModel(j))^2 + (Zdata(i)-ZModel(j))^2);
        distx = sqrt((Xdata(i)-XModel(j))^2);
        disty = sqrt((Ydata(i)-YModel(j))^2);
        distz = sqrt((Zdata(i)-ZModel(j))^2);
        
        if(dist < mindist)
            saveindex = j;
            mindist = dist;
        end
        
        if(distx < mindistx)
            saveindexx = j;
            mindistx = distx;
        end
        
        if(disty < mindisty)
            saveindexy = j;
            mindisty = disty;
        end
        
        if(distz < mindistz)
            saveindexz = j;
            mindistz = distz;
        end
    end
    
    if(mindist < 8)
       percentcnt = percentcnt + 1; 
    end
    
    Xsave(i) = XModel(saveindex);
    Ysave(i) = YModel(saveindex);
    Zsave(i) = ZModel(saveindex);
    
    error(i) = mindist;
    errorx(i) = mindistx;
    errory(i) = mindisty;
    errorz(i) = mindistz;
    
    saveindex = 1;
    mindist = 1000;
    mindistx = 1000;
    mindisty = 1000;
    mindistz = 1000;
end

percentcnt / lengthdata

avgerr = mean(error)

ninefifth = prctile(error,95)
 
percent8 = (nnz(error(:) < 8)/length(error))*100
percent6 = (nnz(error(:) < 6)/length(error))*100
percent4 = (nnz(error(:) < 4)/length(error))*100
percent2 = (nnz(error(:) < 2)/length(error))*100

avgerrx = mean(errorx)

ninefifthx = prctile(errorx,95)
 
percent8x = (nnz(errorx(:) < 8)/length(errorx))*100
percent6x = (nnz(errorx(:) < 6)/length(errorx))*100
percent4x = (nnz(errorx(:) < 4)/length(errorx))*100
percent2x = (nnz(errorx(:) < 2)/length(errorx))*100

avgerry = mean(errory)

ninefifthy = prctile(errory,95)
 
percent8y = (nnz(errory(:) < 8)/length(errory))*100
percent6y = (nnz(errory(:) < 6)/length(errory))*100
percent4y = (nnz(errory(:) < 4)/length(errory))*100
percent2y = (nnz(errory(:) < 2)/length(errory))*100

avgerrz = mean(errorz)

ninefifthz = prctile(errorz,95)
 
percent8z = (nnz(errorz(:) < 8)/length(errorz))*100
percent6z = (nnz(errorz(:) < 6)/length(errorz))*100
percent4z = (nnz(errorz(:) < 4)/length(errorz))*100
percent2z = (nnz(errorz(:) < 2)/length(errorz))*100

% avgerr =  6.9019 For middle arc

%%

avgerr = mean(error)

max(error)

max(Xdata)
max(Ydata)
max(Zdata)
max(X2)
max(Y2)
max(Z2)


min(Xdata)
min(Ydata)
min(Zdata)
min(X2)
min(Y2)
min(Z2)
%avgerr for small arc 4.4215

%avgerr for all 3 = 21.2 mm

