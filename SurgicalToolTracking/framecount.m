%% Read in a bunch of files from a given directory

% target dir where data lives (assume only data files have extension myExt)
dataDir = [pwd '/myCVdataFolder/']
myExt = '.txt';  % assume data files have this extension

% array of fileInof structs
filesInfo = dir(dataDir); 

%remove directories
filesInfo = filesInfo(~[filesInfo.isdir])


fileList = {}; 

%use one file
fileList{1} = 'ToolTrackFrameInfo.txt'

%%

data = load([dataDir fileList{1}]);

[lengther,col] = size(data)

%compute overall metrics
fprintf('Overall Means \n');

total0 = nnz(data(:,2)==0)
total1 = nnz(data(:,2)==1)
total2 = nnz(data(:,2)==2)


percent0 = (total0/lengther)*100
percent1 = (total1/lengther)*100 + (total2/lengther)*100
percent2 = (total2/lengther)*100

% results 
% percent0 =
% 
%    11.9898
% 
% 
% percent1 =
% 
%    88.0102
% 
% 
% percent2 =
% 
%    51.1480