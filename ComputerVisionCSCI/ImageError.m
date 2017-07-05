%% compare two images

ScannerResolution = 23.62; %pixel per mm (600 pixels per inch)
%first crop data image to square 2362x2362 pixels

%load
imtemplate = imread('Minnesota_M.png');
%imdata = imread('Scan_M_Static_3_cropped.png'); %for stationary hand
%imdata = imread('ScanM_Hand_Crop.png'); % for moving model 2D
imdata = imread('Scan_M_Hand_Rod_Pealed2_Crop.png'); % for moving rod 2D

%convert to gray
template = rgb2gray(imtemplate);
data = rgb2gray(imdata);

%show
figure, imshow(template)
figure, imshow(data)

%% first lets upsample the template to match the scan

mmsize = 100;
template = imresize(template,ScannerResolution);

figure, imshow(template)

%%

%threshold
T1 = 125;
BWtemplate = ThresholdImage(template,T1,1);
T2 = 190;
BWdata = ThresholdImage(data,T2,1);

%show
figure, imshow(BWtemplate)
figure, imshow(BWdata)


%%

st = regionprops(BWtemplate, 'Centroid',...
    'MajorAxisLength','MinorAxisLength');
meantemplate = [1,1];
best = 0;
for k = 1:numel(st)
    if(st(k).MajorAxisLength > best)
        best = st(k).MajorAxisLength;
        meantemplate = st(k).Centroid;
    end
end

sd = regionprops(BWdata, 'Centroid',...
    'MajorAxisLength','MinorAxisLength');
meandata = [1,1];
best = 0;
for k = 1:numel(sd)
    if(sd(k).MajorAxisLength > best)
        best = sd(k).MajorAxisLength;
        meandata = sd(k).Centroid;
    end
end
meantemplate
meandata

%center the image
imcenter = size(BWdata)./2;
shifttmeplate = fliplr(round(imcenter - meantemplate)); %(row column)
centertemplate = circshift(BWtemplate,shifttmeplate);
shiftdata =  fliplr(round(imcenter - meandata));%(row column)
centerdata = circshift(BWdata,shiftdata);

%show
figure, imshow(centertemplate)
figure, imshow(centerdata)

%% now we fine tune our matching (keep running this cell till we're happy)

datasize = size(BWdata);
center = round(datasize./2);
shiftdist = [5,5]; %round(datasize ./ 20);
step=1;

theta = -5:1:5;
shiftr = [-shiftdist(1):step:shiftdist(1)];
shiftc = [-shiftdist(2):step:shiftdist(2)];

besterr = mean(mean(abs(double(centerdata - centertemplate))))
bestData = centerdata;
bestParams = [0,0,0];
for th=theta
    th
   for rr = shiftr
       rr
        for cc = shiftc
            rotatedata = imrotate(centerdata,th,'crop');
            shiftdata = circshift(rotatedata,[cc,rr]);
            err = mean(mean(abs(double(shiftdata - centertemplate))));
            if err < besterr
                besterr = err;
                bestData = shiftdata;
                bestParams = [th,cc,rr];
            end
        end
   end
end

bestParams
centerdata = bestData;
figure, imshow(bestData)

%% average total error

differenceim = abs(double(centerdata - centertemplate));

%sum all the wrong locations
errorpixels = differenceim ~= 0;
totalerrorpixels = sum(sum(errorpixels))
totalpixels = size(differenceim,1)*size(differenceim,2)

%average error
averageerror = totalerrorpixels / totalpixels

figure, imshow(differenceim)

%%  Now subtract them

[RR,CC] = size(centertemplate);

diffImage = zeros(RR,CC,3);

truepositive = 0; %inked where we wanted to
falsepositive = 0; %inked where we shouldnt have
truenegative = 0; %didnt ink where we should have
counttrue = 0;
countfalse = 0;
for rr = 1:RR
   for cc = 1:CC
      templatepx =  centertemplate(rr,cc);
      datapx =  centerdata(rr,cc);
      diff = templatepx - datapx;
      if(templatepx >= 255)
         if(diff == 0)
             %inked where we wanted to
             truepositive = truepositive + 1;
             diffImage(rr,cc,2) = 255; %green
         else
             %didnt ink where we should have
             truenegative = truenegative + 1;
             diffImage(rr,cc,3) = 255; %blue
         end
         counttrue = counttrue + 1; %total inks desired
      else
          if(datapx >= 255)
              %inked where we shouldnt have
              falsepositive = falsepositive + 1;
              diffImage(rr,cc,1) = 255; %red
          end
      end
      
      if(datapx >= 255)
          %total inks we deposited
          countfalse = countfalse + 1;
      end
   end
end

truepositive 
falsepositive 
truenegative 
meantruepositive = truepositive / counttrue
meantfalsepositive = falsepositive / countfalse
meantruenegative= truenegative / counttrue

figure, imshow(diffImage)
hold on
h = zeros(3, 1);
h(1) = plot(NaN,NaN,'sg');
hold on
h(2) = plot(NaN,NaN,'sb');
hold on
h(3) = plot(NaN,NaN,'sr');
hold off
legendHandler = legend(h, 'TruePositive','TrueNegative','FalsePositive');
set(legendHandler,'FontSize',32);






