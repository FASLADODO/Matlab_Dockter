%%function [] = ColorIndex(imstring)
%% 
%%Rod Dockter
%%Color Index Work
%%input and preprocessing

clear all;
close all;
clc;

%%getting video
imaqreset;
vid = videoinput('winvideo', 1);
set(vid, 'FramesPerTrigger', Inf);
set(vid, 'ReturnedColorspace', 'rgb')
vid.FrameGrabInterval = 5;

%start video aquisition 
start(vid)
pre = getsnapshot(vid);
[row,col]=size(rgb2gray(pre));
colorimage = pre;
figure('name','initial image')
imshow(colorimage)
%%Declaring variables for histogram work later on
hbin = 12;
sbin = 12;
ibin = 12;
%%using linspace to get nearest bins for each color
hi = linspace(0,1,hbin);
si = linspace(0,1,sbin);
ii = linspace(0,1,ibin);
%%Creating Disk of radius r
radius = (row+col)/7;
D_r = zeros(row,col);
for i = 1:row
   for j = 1:col
      metric = sqrt((i-(row/2))^2+(j-(col/2))^2);
      if metric <radius
         D_r(i,j) = 1; 
      end
   end
end

mainim=colorimage;

%%getting the 'model' 3d color histogram
%%getting the main object matrix which I will assume to be in the center of
%%the image. This will allow me to get the primary color we wish to track.
modelim =rgb2hsi(mainim);
%%for now
rowrange = round((row/2)):1:round((row/2)+(row/5));
colrange = round((col/2)-(col/6)):1:round((col/2)+(col/6));
centerimage = modelim(rowrange,colrange,:);
hcenter = centerimage(:, :, 1);
scenter = centerimage(:, :, 2);
icenter = centerimage(:, :, 3);
htpm = interp1(hi,1:numel(hi),double(hcenter(:)),'nearest');
stpm = interp1(si,1:numel(si),double(scenter(:)),'nearest');
itpm = interp1(ii,1:numel(ii),double(icenter(:)),'nearest');
m3dhist = accumarray([htpm,stpm,itpm],1,[hbin,sbin,ibin]);

%%Getting the color histograms for the image as a whole
%%getting the individual color matrices for image
imagew = rgb2hsi(mainim);
hmat = imagew(:, :, 1);
smat = imagew(:, :, 2);
imat = imagew(:, :, 3);
htpi = interp1(hi,1:numel(hi),double(hmat(:)),'nearest');
stpi = interp1(si,1:numel(si),double(smat(:)),'nearest');
itpi = interp1(ii,1:numel(ii),double(imat(:)),'nearest');
i3dhist = accumarray([htpi,stpi,itpi],1,[hbin,sbin,ibin]);

% Set a loop that stop after 100 frames of aquisition
while(vid.FramesAcquired<=200)
    
    % Get the snapshot of the current frame
    data = getsnapshot(vid);
    
    %%getting ratio histogram (for intersection)
    Rihist = min(m3dhist./i3dhist,1);
    %%getting the histogram color mapping
    workim=rgb2hsi(data);
    him = workim(:, :, 1);
    sim = workim(:, :, 2);
    iim = workim(:, :, 3);
    hind = interp1(hi,1:numel(hi),double(him),'nearest');
    sind = interp1(si,1:numel(si),double(sim),'nearest');
    iind = interp1(ii,1:numel(ii),double(iim),'nearest');
    bidx = sub2ind(size(Rihist),hind,sind,iind);
    b_r=Rihist(bidx);
    b_c = filter2(D_r,b_r);
    [num idx] = max(b_c(:));
    [xc,yc] = ind2sub(size(b_c),idx);
    disp(xc)
    disp(yc)
    disp(num)
    metnum = max(b_c(:));
    area = sum(b_c(:)>metnum/2.3);
    R = sqrt(area/pi);
    %This is a rough estimate of the radius 
    x=0:0.01:1; %degree 
    hold on;
    plot(yc+R*cos(2*pi*x),xc+R*sin(2*pi*x), 'w')
    hold off;

end



stop(vid);

flushdata(vid);

%%
%%end