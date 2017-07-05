
% Read in image and display color scale
inputpic = imread('Monkey.JPG');

%get image size and channels
[row,col,channels]=size(inputpic);


figure(1)
imshow(inputpic)
title('Orginal Image')


%% matlab built in grayscale

gray = rgb2gray(inputpic);
figure(2)
imshow(gray)
title('Automatic grayscale')

minval = min(min(gray))
maxval = max(max(gray))

%% extract each channel

redchannel = inputpic(:,:,1);
greenchannel = inputpic(:,:,2);
bluechannel = inputpic(:,:,3);


figure(3)
imshow(bluechannel)
title('Individual blue channel')


%% manual grayscale

%create single channel matrix of zeros
mangray = zeros(row,col);
%convert back into unsigned integer
mangray = uint8(mangray);

for i=1:row
    for j=1:col
        %do math in double (uint8 is limited to 255 max value)
        pixelval = (double(redchannel(i,j)) + double(greenchannel(i,j)) + double(bluechannel(i,j))) / 3.0 ;
        mangray(i,j) = pixelval;
    end
end



figure(4)
imshow(mangray)
title('Manual grayscale')

%% Threshold

%create single channel matrix of zeros
threshold = zeros(row,col);
threshold = uint8(threshold);
threshval = 120;

for i=1:row
    for j=1:col
        pixelval = mangray(i,j);
        if(pixelval > threshval)
           threshold(i,j) = 255; 
        end
    end
end

figure(5)
imshow(threshold)
title('Manual Thresholding')


%% Flip Colors

%create three channel matrix of zeros
Flipped = zeros(row,col,3);
Flipped = uint8(Flipped); %convert

for i=1:row
    for j=1:col
        %map rgb to gbr
        Flipped(i,j,2) = redchannel(i,j);
        Flipped(i,j,3) = greenchannel(i,j);
        Flipped(i,j,1) = bluechannel(i,j);
    end
end


figure(6)
imshow(Flipped)
title('RGB -> GBR')

%% Invert Colors

%create three channel matrix of zeros
Inverted = zeros(row,col,3);
Inverted = uint8(Inverted); %convert

for i=1:row
    for j=1:col
        %invert magnitude of each channel
        Inverted(i,j,1) = 255 - redchannel(i,j);
        Inverted(i,j,2) = 255 - greenchannel(i,j);
        Inverted(i,j,3) = 255 - bluechannel(i,j);
    end
end


figure(7)
imshow(Inverted)
title('Inverted Color')


%% Edge Detection


edgepic = zeros(row,col);
edgepic = uint8(edgepic);

picg = threshold; %gray or individual channel

for i=2:row-1
    for j=2:col-1
        %%Produce gradient terms with sobel and roberts operator
        deltas1=double( picg(i-1,j+1)+2*picg(i,j+1)+picg(i+1,j+1)-picg(i-1,j-1)-2*picg(i,j-1)-picg(i+1,j-1) );
        deltas2= double( picg(i-1,j-1)+2*picg(i-1,j)+picg(i-1,j+1)-picg(i+1,j-1)-2*picg(i+1,j)-picg(i+1,j+1) );
        mags=sqrt((deltas1^2)+(deltas2^2));
        if mags > 25
            edgepic(i,j)=255;
        end
    end
end

figure(8)
imshow(edgepic)
title('Edge Detection')
