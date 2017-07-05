
dataDir = [pwd '/myCVdataFolder/']
myExt = '.txt';  

filesInfo = dir(dataDir); 

filesInfo = filesInfo(~[filesInfo.isdir])

fileList = 'DepthCalibration.txt'

data = load([dataDir fileList]);
size(data)


%%
% Get Middle Point For X and Y from Each Channel
X = (data(:,1) + data(:,3)) / 2.0;
Y = (data(:,2) + data(:,4)) / 2.0;

% Disparity
disp = data(:,5);

datad = [X,Y,disp];

realX = data(:,6);
realY = data(:,7);
realZ = data(:,8);

modeld = [realX,realY,realZ];

[R,T] = icp(modeld,datad);
