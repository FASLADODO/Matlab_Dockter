%%function [dlist,thetalist] = HoughLineFit(imagestring,votethreshold,dbinsize,thetabinsize)
%% 
%%Rod Dockter
%%Programming Assignment 2
%%Csci 5561 Computer Vision
%%Again done in cell mode but can easily be function ready

%%The input arguments consist of:
%%the name of the image in question
%%the threshold for edge detection in the image
%%the minimum threshold for vote counting
%%the size of each bin in the accumulator array for D
%%the size of each bin in the accumulator array for theta
imagestring = 'line1.jpg';
votethreshold= 100;
dbinsize = 5;
thetabinsize = 1;


%%This function will perform the hough transform for line fitting to the
%%line1 image provided for this assignment. It will do this via the classic
%%hough transform method of building up the accumulator array
%%of all possible lines and then counting votes for lines.
%%I will use the thinning algorithm from programming assignment
%%so that the lines the accumulation deals with are only 1 pixel thick.

%%Calling the edgedetect/thinning algorithm
[thinnedpic] = EdgeDetectDockter(imagestring,0.38,4);

%%Getting size of image
[row,col]=size(thinnedpic);

%%Now for the accumulator array
%%D will be rows of array
%%theta will be columns of array
%%I will define the asix to just follow the axis of the pixels in the image
%%so that plotting will be easier later on
%%The range of the d will therefore need to go from the bottom left corner
%%to the top right corner meaning the range will be  +/- sqrt(row^2+ col^2)
%%and the theta range will be 0 - 180

%%Define range of theta and r accumulator values
dpoints = -sqrt(row^2+col^2):dbinsize:sqrt(row^2+col^2);
thetapoints = 0:thetabinsize:180-thetabinsize;

%%initializing accumulator to zeros
Hacc=zeros(length(dpoints),length(thetapoints));

%%Creating array of x and y coordinates of edge points in image
xcoord = [];
ycoord = [];
index = 1;
for i=1:row
    for j=1:col
        if thinnedpic(i,j) == 1
            xcoord(index) = j;
            ycoord(index) = i;
            index=index+1;
        end
    end
end

%%%%%%%%%%%%%%%%%This section merely populates the the accumulator array
%%%%%%%%%%%%%%%%%based of the equationr =x*cos(theta)+y*sin(theta)

%%Looping through edge points
%%Then looping through theta values
%%Finding D and populating in accumulator array
for k = 1:length(ycoord)
    %%For each new edge point I reset the theta index and loops through my
    %%desired theta values
    indextheta = 1;
    for th=thetapoints
        d = xcoord(k)*cosd(th)+ycoord(k)*sind(th);
        %%Another loop to find out the desired index in the accumulator
        %%array for the D value
        min_d = 10000;
        for h=1:length(dpoints)
            if abs(d-dpoints(h)) < min_d
               indexd = h;
               min_d = abs(d-dpoints(h));
            end
        end
        %%Once proper index value is found, increase accumulator array
        %%value
        Hacc(indexd,indextheta) = Hacc(indexd,indextheta) + 1;
        %%increase theta indexer
        indextheta = indextheta+1;
    end
end

%%%%%%%%%%%The next section goes through the Accumulator Array and finds
%%%%%%%%%%%the largest number of votes.

%%preallocate
dparam = [];
thetaparam=[];
%%looping through and finding the desired points
%%In order to ensure I only use local maxima I will check each neighbor
%%for points above the basic threshold to see if the point is larger than
%%all it's neighbors in the 5x5 square surrounding it.
filltick = 1;
for i=1:length(dpoints)
    for j=1:length(thetapoints)
        if Hacc(i,j) > votethreshold
            neighbori = -2:2;
            neighborj = -2:2;
            neighbors = zeros(1,48);
            ntick=1;
            %%checking all the neighbors (in range) within 5 of the
            %%accumulator vote in question
            %%IF any of the nearby votes are larger I can ignore a point.
            for ki = neighbori
                for kj = neighborj
                    tempi = i+ki;
                    tempj = j+kj;
                    if(tempi ~= 0) && (tempj ~= 0)
                        if (tempi <= length(dpoints)) && (tempj <= length(dpoints))
                            neighbors(ntick) = Hacc(tempi,tempj);
                            ntick=ntick+1;
                        end
                    end
                end
            end
            %%If the accumulator point has enough votes I record the theta
            %%and d values minus half the binsize to center the values in
            %%bin
            %%Checking local maxima-ship (this is a made up word)
            if (any(Hacc(i,j) < neighbors) == 0)
                dparam(filltick) = dpoints(i)-dbinsize/2;
                thetaparam(filltick) = thetapoints(j)-thetabinsize/2;
                filltick=filltick+1;
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%This final section is for return arguments and for the
%%%%%%%%%%%%%%%%%%plotting of figures.
dlist=dparam;
thetalist=thetaparam;

figure('Name','Thinned Image');
imshow(thinnedpic)

%%plotting the lines over the image
%%Plots all the resultant lines using ezplot
figure('Name','Image with Hough Lines');
readin = imread('line1.jpg');
imshow(readin);
hold on
for i = 1:length(dlist)
    fun = @(x,y) x*cosd(thetalist(i))+y*sind(thetalist(i))-dlist(i);
    ezplot(fun,[1 600 1 800])
    hold on
end
hold off


% %%Plotting the accumulator for visualization
figure('Name','Accumulator Array');
mesh(Hacc)
%% 
%%end





