%Take Home Final
%Problem 8
%Rod Dockter
%Grassfire Transform
%Code is similar to what I used for programming assignment 3, with a few
%adjustments
%I just used the ellipse image, I figured that was sufficient for
%demonstration purposes, the image can be changed to any other one.

%Resultant images can be found at

%https://sites.google.com/a/umn.edu/umnworkdockter/grassfiretransform

clear all;
close all;
clc;


%%Reading in image data and converting to black and white
picnamestring = 'ellipse1.jpg';
inputpic = imread(picnamestring);
%Getting the black and white image
inputpicb =im2bw(inputpic);
%putting in double type so range is 0 -1
picd = im2double(inputpicb);
%getting size
[row,col]=size(picd);

%displying original image in bw
figure('Name','input image');
imshow(inputpic)

%preallocating the grassfire matrix
grassfire = zeros(row,col);

%this first loop scans through the image matrix from top left to bottom
%right, if a pixel is black, it's position in the grassfire matrix is then
%given a value equal to the minimum of it's neighbors values to the north
%and west plus one. So that as we sweep down and to the right through and
%object the distance values in the grassfire will get larger and larger.
%If only this loop was used, only the distances to the top and left object
%edge wouold be found.
for i = 2:row-1
    for j = 2:col-1
        val = picd(i,j);
        %if pixel is object pixel
        if val == 0
            %%getting all neighbors (N+W)
            gfneighbor = [grassfire(i-1,j-1),grassfire(i-1,j),grassfire(i-1,j+1),grassfire(i,j-1),];
            %setting the value to the minimum neighbors value + 1
            grassfire(i,j) = min(gfneighbor)+1;
        end
    end
end

%Now doing the second pass the oppositie direction through the image
%here we scan from the bottom right corner to the top left. Using the
%partially formulated grassfire matrix from the first pass we check if a
%pixel is and object pixel and if it we replace the distance value in the
%grassfire matrix with the value of its minimum neighbor but now using the
%neighbors to the right and down (E+S).
%The effect of this is that distance of each object pixel to the right and
%down edges of it's object are also considered.
%Notice, to go through the opposite direction of the grassfire matrix I use
%the decrement operator (-1)
for i = row-1:-1:2
    for j = col-1:-1:2
        val = picd(i,j);
        %check if pixel is object
        if val == 0
            %%getting all neighbors (S+E)
            gfneighbor = [grassfire(i+1,j-1),grassfire(i+1,j),grassfire(i+1,j+1),grassfire(i,j+1)];
            %Setting the grassfire value to the minimum of it's existing
            %value (if the distance was already correct) and the minimum
            %(S+E) neighbor + 1
            grassfire(i,j) = min(grassfire(i,j),min(gfneighbor)+1);
        end
    end
end

%Notice both the above loops had to go from 2 through col/row - 1 since I
%was considering neighbors (just being thorough)

%%displaying grassfire visualization

%this factor helps add definition to the skeleton
sfactor = 1.1; 
%making a skeleton matrix which is just the grassfire matrix divided by the
%max distance (center of largest object in grassfire.
%This scales the gradient so that it is between 0-1 (white and black)
%corresponding to how far the pixel is from the edge.
skeleton = grassfire/(sfactor*max(grassfire(:)));

%Displaying the skeleton image
figure('Name','grassfire transform')
imshow(skeleton)



%%