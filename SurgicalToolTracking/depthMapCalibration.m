
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
X = (data(:,1) + data(:,3)) / 2.0;
Y = (data(:,2) + data(:,4)) / 2.0;

disp = data(:,5);

indata = [disp,X-xcenter,ycenter-Y];

ZData = data(:,8);

%%

guessA = [31992,-1.117,0.0081,0.0256,0.3];

realA = nlinfit(indata,ZData,@depthCalPower,guessA);

fprintf('%.8f ', realA)

%%
% recent
% realA =
% 
%   712.8527   -0.0137   -0.0163   -0.0708    7.1985
% for Z = A*e^B*disp + C*(X-Xc) + D*(Y-Yc) + E

% realA =
% 
%   709.1749   -0.0132   -0.0163   -0.0707
% Z = A*e^B*disp + C*(X-Xc) + D*(Y-Yc) + E
%%

indataX = [data(:,8),X-xcenter];

XData = data(:,6);

guessB = [0.0016,0.3643];

realB = nlinfit(indataX,XData,@XYCal,guessB);

% fprintf('%.8f ', realB)

realB

%%

indataY = [data(:,8),ycenter-Y];

YData = data(:,7);

guessC = [0.0016,0.0014];

realC = nlinfit(indataY,YData,@XYCal,guessC);

% fprintf('%.8f ', realC)

realC

