%%function [] = ColorIndex(imstring)
%% 
%%Rod Dockter
%%Color Index for static image

clear all;
close all;
clc;

%%Reading in image data and converting to greyscale
% picnamestring1 = 'cluttered.jpg';
% 
% inputpic1 = imread(picnamestring1);
% 
% colorimage1 = im2double(inputpic1);
% [row1,col1]=size(rgb2gray(colorimage1));
% 
% figure('name','lighter image')
% imshow(colorimage1)


%%Reading in image data and converting to greyscale
picnamestring = 'cluttered.jpg';

inputpic = imread(picnamestring);

colorimage = im2double(inputpic);
[row,col]=size(rgb2gray(colorimage));

figure('name','cluttered darker image')
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
radius = (row+col)/12;
D_r = zeros(row,col);
for i = 1:row
   for j = 1:col
      metric = sqrt((i-(row/2))^2+(j-(col/2))^2);
      if metric <radius
         D_r(i,j) = 1; 
      end
   end
end

mainim1=colorimage1;
mainim = colorimage;

time1 = cputime;
%%getting the 'model' 3d color histogram
%%getting the main object matrix which I will assume to be in the center of
%%the image. This will allow me to get the primary color we wish to track.
modelim =rgb2hsi(mainim1);
%%for now
modelcx = 262;
modelcy = 150;
modelr = 10;
rowrange = round(modelcy-modelr):1:round(modelcy+modelr);
colrange = round(modelcx-modelr):1:round(modelcx+modelr);
centerimage = modelim(rowrange,colrange,:);
hcenter = centerimage(:, :, 1);
scenter = centerimage(:, :, 2);
icenter = centerimage(:, :, 3);
htpm = interp1(hi,1:numel(hi),double(hcenter(:)),'nearest');
stpm = interp1(si,1:numel(si),double(scenter(:)),'nearest');
itpm = interp1(ii,1:numel(ii),double(icenter(:)),'nearest');
m3dhist = accumarray([htpm,stpm,itpm],1,[hbin,sbin,ibin]);


comptime = cputime - time1;

disp(comptime)

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

    
% Get the snapshot of the current frame
data = mainim;
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
b_c = real(fftshift(ifft2(fft2(b_r).*fft2(D_r))));
[num idx] = max(b_c(:));
[rowc,colc] = ind2sub(size(b_c),idx);
disp(rowc)
disp(colc)
R = 20;
%This is a rough estimate of the radius 
hold on;
%now adding circle depicting backprojected location
rectangle('Position',[colc-R/2,rowc-R/2,R,R],'EdgeColor','red')
hold off




%%
%%end