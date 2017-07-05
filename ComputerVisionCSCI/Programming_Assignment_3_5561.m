%%
%%Programming Assignment # 3
%%Image Recognition
%%Rod Dockter
%%Csci 5561
%%Notice this done in cell mode (could easily be made into full function)

%%The purpose of this code is to perform (in logical chunks) a histogram
%%threshold determination of an image filled with different objects, then
%%to create an image of object regions so that we can perform connected
%%region extraction on the different unique regions
%%Then finally we will compute certain 'blob' statistics for each region
%%In addition to the required work I added a grassfire transform for
%%potential future use and a run compression of the binary image.

%%Notice that due to shadows, some of my objects get connected. This is
%%less than ideal, but every thing else works fine.

clear all;
close all;
clc;


%%Reading in image data and converting to greyscale
picnamestring = 'reg4.jpg';
disp('image in use: ')
picnamestring
inputpic = imread(picnamestring);
%%uint8 version
picg = im2uint8(inputpic);
%%double version
picd = im2double(inputpic);
[row,col]=size(picg);
figure('Name','input image');
imshow(picd)
%%Preallocating histogram space
%%Note: hisotogram will be 1- 256 even thought intensity is 0 -255 for
%%indexing reasons
inthist = zeros(1,256);

%%
%%filling hisotgram with intensity values
for i = 1:row
   for j = 1:col
      %%Have to increment each bin by 1 since intensity is zero based
      intensity = double(round(picg(i,j) + 1));
      inthist(1,intensity) = inthist(intensity) + 1;
   end
end
figure('Name','normal histogram');
plot(inthist)

%%
%%Smoothing the histogram
smhist = zeros(1,256);
%%I do this by summing each bin with its 4 closest neighbors
for i = 1:256
    %%Have to include cases for ends of histogram
    if i==1
       smhist(i) = inthist(i)+inthist(i+1);  
    elseif i == 2
        smhist(i) = inthist(i-1)+inthist(i)+inthist(i+1);
    elseif i == 255
        smhist(i)= inthist(i-1)+inthist(i)+inthist(i+1);
    elseif i == 256
        smhist(i) = inthist(i)+inthist(i-1); 
    else
        %%normal smoothing portion
        smhist(i) = inthist(i-2)+inthist(i-1)+inthist(i)+inthist(i+1)+inthist(i+2);
    end
end
figure('Name','smoothed histogram');
plot(smhist)

%%
%%getting the threshold using otsu method
%%Notice this code was already written for my edge detection program
%%I just reused it.

%%Initialize various values for computation of these
totalpixel = row*col;
BCVmax=0;
magthreshold=0;
directionr=0;
directions=0;
Wbnum=0;
Wbden=totalpixel;
Wb=0;
Wfnum=0;
Wfden=totalpixel;
Wf=0;
Ubnum=0;
Ubden=0;
Ub=0;
Ufnum=0;
Ufden=0;
Uf=0;
for k=1:256
    weightedtotal=(k-1)*(smhist(k)-1);%%has to be decreased by 1 from hist
end
%%Loop for threshold determination (Maximizing between class variance)
for m=1:1:256
    %%Determining the several values for between class variances
    Wbnum=Wbnum+(smhist(m)-1);
    Wfnum=totalpixel-Wbnum;
    Wb=Wbnum/Wbden;
    Wf=Wfnum/Wfden;
    if (Wfnum ~= 0) && (Wbnum~=0)
        Ubnum=Ubnum + m*(smhist(m)-1);
        Ubden=Wbnum;
        Ub=Ubnum/Ubden;
        Ufnum=weightedtotal-Ubnum;
        Ufden=Wfnum;
        Uf=Ufnum/Ufden;
    end
    %%preventing divide by zero
    if (Wfnum == 0) || (Wbnum == 0)
        Ub=0;
        Uf=0;
    end
    %%Between Class Variance
    BCV=Wb*Wf*((Ub-Uf)^2);
    %%Maximize
    if BCV>BCVmax
       BCVmax=BCV;
       magthreshold=m;
    end
end
%%printing the threshold
%%The actual value I get ended up including shadows, so I lowered it by a
%%small percent
magthreshold = magthreshold * 0.9;
%%converting to a corresponding threshold for the double 0-1 scale
doublethresh = magthreshold/255;

%%
%%seperating objects from background
%%this would be inverted if foreground was white objects
objectmat = zeros(row,col);
objectim = mat2gray(uint8(objectmat));
for i = 1:row
   for j = 1:col
       val = picg(i,j);
       %%If object is dark enough it becomes black
       if val <= magthreshold
           objectim(i,j) = 0;
       else
          %%If not you make it white
          objectim(i,j) = 255; 
       end
   end
end
figure('Name','Objects vs background image');
imshow(objectim)

%%
%%Determing if black or white object are the target objects
%%assuming there will be less pixels of the color corresponding to the
%%target object
blacktally = 0;
whitetally = 0;
for i = 1:row
   for j =1:col
       val = objectim(i,j);
       %%Incrementing one or the other 
       if val == 0
          blacktally = blacktally + 1; 
       end
       if val == 255
          whitetally = whitetally + 1; 
       end
   end
end
%%determining which objects there are fewer of, setting desired objects
if blacktally > whitetally
   obval = 255; 
end
if blacktally <= whitetally
   obval = 0; 
end
disp('objects will have value: ')
obval

%%
%%Performing the  connected region extraction using the raster scan
%%I use the 8 connectivity for this
%%Note I am ignoring any objects on the border and focusing on the more
%%central part of the image (since I know this is where the objects are
regionmat = zeros(row,col);
regionindex = 1; %%will get incremented to start new partial region

for i = 5:row-5
    for j = 5:col-20
        val = objectim(i,j);
        %%checking if object is of the b/w we want
        if val == obval
            neighbors = [regionmat(i,j-1),regionmat(i-1,j-1),regionmat(i-1,j),regionmat(i-1,j+1)];
            %%checking all previous region labels
            if sum(neighbors) == 0
                %%if no previous labels (from N/W) make new label
                regionmat(i,j) = regionindex;
                regionindex = regionindex + 1;
            end
            if sum(neighbors) ~= 0
                propind=min(neighbors(neighbors~=0));
                %%propogating previous lables (if present) using minimum
                %%(nonzero) neighbor's label
                regionmat(i,j) = propind;
            end
        end
    end
end
%%resultant index will be the highest label value
%%which we will need
maxregion = regionindex -1;

%%
%%Now getting the equivalence chart for merging the regions
%%Top row will be index, second row will be lowest label in region
eqchart = 1:maxregion;
for i = 2:row-1
    for j = 2:col-1
        val = regionmat(i,j);
        %%if the label is not a zero (background)
        if val ~= 0
            %%finding the labels of all the neighbors
            labneighbor = [regionmat(i-1,j-1),regionmat(i-1,j),regionmat(i-1,j+1),regionmat(i,j-1),regionmat(i,j+1),regionmat(i+1,j-1),regionmat(i+1,j),regionmat(i+1,j+1)];
            minnz = min(labneighbor(labneighbor~=0)); %%finding the minimum non zero neighbor
            eqchart(val) = min([val,eqchart(minnz),eqchart(eqchart(val))]);
            %%this sets the equivilant chart value for the current label
            %%to the minimum of its actual label number, its smallest
            %%neighbor, or its smallest neighbor's equivilant chart value
        end
    end
end

%%Reordering equivalence chart so labels will go from 1 -> total number of
%%regions
%%To do this I utilize the unique() function in matlab which returns an
%%array of the unique labels in the equivalence chart array then use the
%%find() indices function to determine the indices of the unique values
%%corresponding to each value in the equivalence chart.
%%This results in the values in my new equivalence chart ranging from 1 to
%%the number of total different regions (this is more useful than the
%%leftover large label values from the first scan)
totallabels = unique(eqchart);
eqlabel = zeros(1,maxregion);
for i = 1:maxregion
    %%getting new equivalence charts in order
    eqlabel(i) = find(totallabels==eqchart(i));
end

%%so the total number of regions I have is:
disp('the total number of different regions is')
totalregions = max(eqlabel)

%%
%%Now I use the equivalence chart to take a second scan of the object image
%%and merge all the regions to lowest region label
regionlabel = zeros(row,col);
for i = 2:row-1
    for j = 2:col-1
        val = regionmat(i,j);
        if val ~= 0
            regionlabel(i,j)=eqlabel(val);
        end
    end
end
%%
%%Various plotting work
%%Diplaying the different regions as different colors
colorlabel = zeros(row,col,3, 'uint8');
%%Color maps didnt seem to work for an rgb matrix, so I made up colors,
%%this is why the colors end up looking ugly
rcolors = randi([50,230],1,totalregions);
bcolors = randi([50,230],1,totalregions);
gcolors = randi([50,230],1,totalregions);
for i = 1:row
    for j = 1:col
        val = regionlabel(i,j); 
        if val ~= 0
           colorlabel(i,j,1) = rcolors(val);
           colorlabel(i,j,2) = bcolors(val); 
           colorlabel(i,j,3) = gcolors(val); 
        else
           colorlabel(i,j,1) = 255;
           colorlabel(i,j,2) = 255; 
           colorlabel(i,j,3) = 255; 
        end
    end
end
figure('Name','Color Coded Regions');
imshow(colorlabel)

%%making an individual matrix 
indivmats = zeros(row,col,totalregions);
for k = 1:totalregions
    for i = 1:row
        for j = 1:col
            val = regionlabel(i,j);
            if val == k
               indivmats(i,j,k) = 0; 
            else
               indivmats(i,j,k) = 255;
            end
        end
    end
end

plotlength = round(totalregions/2);
figure('Name','Individual regions');
for i = 1:totalregions
    subplot(plotlength,2,i);
        subimage(indivmats(:,:,i))
        title(['Plot of Region # ',num2str(i)])
end

%%
%%The next sections are for getting blob statistics
%%Notice whenever I print statistics, they will be in order from
%%1->totalregions as defined by the region numbers in the above graphs

%%getting area of each region in area
Area = zeros(1,totalregions);
for k = 1:totalregions
    %%simply finding the sum of the times each label occurs in the label
    %%matrix
    Area(k)=sum(regionlabel(:)==k);
end
disp('the areas for each region are')
Area

%%
%%Getting minimum bounding rectangle (MBR)
MBR = zeros(4,totalregions);
%%data will be stored in 4xn array
%%1st row = minX, 2nd row = minY, 3rd row = maxX, 4th row = maxY
for k = 1:totalregions
    %%using temp values to find minimum for each label
    mintx = col;
    minty = row;
    maxtx = 0;
    maxty = 0;
    minycorx = 0;
    for i = 1:row
        for j = 1:col
            val = regionlabel(i,j);
            if val == k
                %%updating minimum/maximum values for MBR
                if j > maxtx
                   maxtx = j;
                end
                if j < mintx
                   mintx = j;
                end
                if i > maxty
                   maxty = i;
                end
                if i < minty
                   minty = i;
                   minycorx = j;
                end
            end
        end
    end
    %%storing resultant data in array
    MBR(1,k) = mintx;
    MBR(2,k) = minty;
    MBR(3,k) = maxtx;
    MBR(4,k) = maxty;
end

%%Displaying arrays of MBR data
disp('the minimum X coordinates of the MBRs are')
xminMBR = MBR(1,:)
disp('the minimum Y coordinates of the MBRs are')
yminMBR = MBR(2,:)
disp('the maximum X coordinates of the MBRs are')
xmaxMBR = MBR(3,:)
disp('the maximumY coordinates of the MBRs are')
ymaxMBR = MBR(4,:)

%%Plotting the square around the original image for a better visualization
%%I use a color map to have a random different color for each MBR
cmap = hsv(totalregions);
figure('Name','MBR around image')
imshow(picg)
hold on
for i = 1:totalregions
    %%Plotting each square (one corner is repeated in order to complete
    %%square)
    plot([MBR(1,i),MBR(1,i),MBR(3,i),MBR(3,i),MBR(1,i)],[MBR(2,i),MBR(4,i),MBR(4,i),MBR(2,i),MBR(2,i)],'-s','Color',cmap(i,:))
    hold on
end
hold off

%%
%%Getting centroid data
%%all coordinate information will have the origin at the top left corner
xtally = zeros(1,totalregions);
ytally = zeros(1,totalregions);

for k = 1:totalregions
    for i = 1:row
        for j = 1:col
            val = regionlabel(i,j);
            if val == k
                xtally(k) = xtally(k) + j; 
                ytally(k) = ytally(k) + i;
            end
        end
    end
end
xcentroid = xtally ./ Area;
ycentroid = ytally ./ Area;

%%Plotting the centroids over the original image for visualization
%%used a colormap again to have a random value
cmap = hsv(totalregions);
figure('Name','Centroids on image')
imshow(picg)
hold on
for i = 1:totalregions
    %%Plotting each square (one corner is repeated in order to complete
    %%square)
    plot(xcentroid(i),ycentroid(i),'s','Color',cmap(i,:))
    hold on
end
hold off
disp('the x centroid coordinates are')
xcentroid
disp('the y centroid coordinates are')
ycentroid


%%
%%Finding the holes using an psuedo reverse region labeling
%%I just took the region labeling code and augmented it

regionholemat = zeros(row,col);
regionhindex = 1; %%will get incremented to start new partial region
for i = 2:row-1
    for j = 2:col-1
        hval = objectim(i,j);
        %%checking if pixel is a hole (white is this case)
        if hval == 255
            hneighbors = [regionholemat(i,j-1),regionholemat(i-1,j-1),regionholemat(i-1,j),regionholemat(i-1,j+1)];
            %%checking all previous region labels
            if sum(hneighbors) == 0
                %%if no previous labels (from N/W) make new label
                regionholemat(i,j) = regionhindex;
                regionhindex = regionhindex + 1;
            end
            if sum(hneighbors) ~= 0
                hpropind=min(hneighbors(hneighbors~=0));
                %%propogating previous lables (if present) using minimum
                %%(nonzero) neighbor's label
                regionholemat(i,j) = hpropind;
            end
        end
    end
end
%%resultant index 
maxhregion = regionhindex -1;
heqchart = 1:maxhregion;
for i = 2:row-1
    for j = 2:col-1
        hval = regionholemat(i,j);
        %%if the label is not a zero (objects in this case)
        if hval ~= 0
            %%finding the labels of all the neighbors
            labhneighbor = [regionholemat(i-1,j-1),regionholemat(i-1,j),regionholemat(i-1,j+1),regionholemat(i,j-1),regionholemat(i,j+1),regionholemat(i+1,j-1),regionholemat(i+1,j),regionholemat(i+1,j+1)];
            minnz = min(labhneighbor(labhneighbor~=0)); %%finding the minimum non zero neighbor
            heqchart(hval) = min([hval,heqchart(minnz),heqchart(heqchart(hval))]);
        end
    end
end
%%Reordering hole equivalence chart 
totalhlabels = unique(heqchart);
heqlabel = zeros(1,maxhregion);
for i = 1:maxhregion
    %%getting new equivalence charts in order
    heqlabel(i) = find(totalhlabels==heqchart(i));
end

%%Now I use the equivalence chart to take a second scan of the object image
%%and merge all the connected holes
regionholelabel = zeros(row,col);
for i = 2:row-1
    for j = 2:col-1
        hval = regionholemat(i,j);
        if hval ~= 0
            regionholelabel(i,j)=heqlabel(hval);
        end
    end
end

%%Now I find any hole region that touches a boundary and remove that as a hole
%%region
for i = 2:row-1
    for j = 2:col-1
        hval = regionholelabel(i,j);
        if hval ~= 0
            %%I leave some room for error in 'touching' boundary of image
            if (i >= row-1) || (i <= 10) || (j >=col-10) || (j <= 10)
                %%if that 'hole' region is near boundary I remove that
                %%label from the whole image
                regionholelabel(regionholelabel == hval)=0;
            end
        end
    end
end

%%Now to count up the number of holes and their areas
%%Again using unique
uniarray = unique(regionholelabel);
uniarray = uniarray(uniarray ~= 0);
%%getting number of unique hole labels (non zero since those arent holes)
numberofholes = length(uniarray);
areaofholes = zeros(1,numberofholes);
%%preallocating
for k = 1:numberofholes
   regval =  uniarray(k);
   if regval ~= 0
       %%for each hole label value, I count up the number of pixels in hole
       %%label image with that value and return that as the area
       areaofholes(k) = sum(regionholelabel(:) == regval);
   end
end
%%printing the results
disp('number of holes is: ')
numberofholes
disp('area of holes are: ')
areaofholes

%%showing labels of holes in different colors
%%This wasnt required so I figured it was okay to use the rgb2label method
%%just for me to visualize
figure('Name','colored unique hole regions');
labelholeim = label2rgb(regionholelabel);
imshow(labelholeim)

%%
%%Getting Perimeter statistics for each region
%%First have to get a thinned image
%%To do this I will write a psuedo edge detection for each 
labeledge = zeros(row,col);
%%I will go through and store any region label that has a neighboring 0 as
%%an edge point. This is not perfect but helps to prevent (further) merging
%%of regions.
%%I also use this loop to find the minimum edge point  for perimeter
%%searching
minpoint = zeros(2,totalregions);
minpoint(:) = max(row,col); %%sets to a max value to find min point
for k = 1:totalregions
    for i = 1:row
        for j = 1:col
            val = regionlabel(i,j);
            if val == k
                %%getting all neighbors
                rneighbor = [regionlabel(i-1,j-1),regionlabel(i-1,j),regionlabel(i-1,j+1),regionlabel(i,j-1),regionlabel(i,j+1),regionlabel(i+1,j-1),regionlabel(i+1,j),regionlabel(i+1,j+1)];
                %%checking if there are any zeros (non regions) in the 8
                %%connectvity of the region point in question
                if any(rneighbor == 0)
                    labeledge(i,j) = k;
                    if i <= minpoint(2,k)
                        %%recording minimum edge point.
                        minpoint(2,k)=i;
                        minpoint(1,k)=j;
                    end
                end
            end
        end
    end
end

%%The above loop also considers inner holes. Therefore I go through a
%%second loop to remove extra holes by checking if they are connected to
%%hole label regions
perimlabel = zeros(row,col);
for k = 1:totalregions
    for i = 2:row-1
        for j = 2:col-1
            val = labeledge(i,j);
            if val == k
                holeneighbors = [regionholelabel(i,j),regionholelabel(i-1,j-1),regionholelabel(i-1,j),regionholelabel(i-1,j+1),regionholelabel(i,j-1),regionholelabel(i,j+1), regionholelabel(i+1,j-1),regionholelabel(i+1,j),regionholelabel(i+1,j+1)];
                if (any(holeneighbors(:)~=0)) 
                    perimlabel(i,j) = 0;
                else
                    perimlabel(i,j) = k;
                end
            end
        end
    end
end

%%This is just a loop to take the 'thinned' region image and turn it into a
%%black and white image of these edges for visualization
edgemat = zeros(row,col);
for i = 1:row
    for j = 1:col
        val = perimlabel(i,j);
        %%turning any label into a black pixels
        if val ~= 0
           edgemat(i,j) = 0; 
        else
           edgemat(i,j) = 255; 
        end
    end
end
figure('Name','perimeter pixels')
imshow(perimlabel)
hold on
for i = 1:totalregions
    %%Plotting each square (one corner is repeated in order to complete
    %%square)
    plot(minpoint(1,i),minpoint(2,i),'s','Color',cmap(i,:))
    hold on
end

%%Then finally I can loop through and find the total number of each region
%%in the thinned edge image
perimeters = zeros(1,totalregions);
for i =1:totalregions
    perimeters(i) = sum(perimlabel(:) == i);
end
disp('the length of perimeter for each region is')
perimeters

%%Therefore we can compute the elongation as P^2/A
psqr = perimeters.*perimeters;
elongation = psqr./Area;
disp('The elongations are ')
elongation



%%
%%EXTRA STUFF
%%turning image into run compression of the image (extra)
%%Initially I had hoped to use this for hole statistics, but that didnt
%%work for me.
prevval = 0;
runmat(1:row,1,1)=0;
for i =1:row
   prevval = objectim(i,1);
   nind = 1;
   for j = 1:col
      currentval = objectim(i,j);
      if currentval == prevval
          runmat(i,nind,1) = runmat(i,nind,1)+1;
          runmat(i,nind,2) = currentval;
      else
          nind = nind+1;
          runmat(i,nind,1) = 1;
          runmat(i,nind,2)=currentval;
          prevval = objectim(i,j);
      end
   end
end

%%
%%For fun I implemented a grassfire transform and plotted the results
%%using grassfire transform
gftran = zeros(row,col,totalregions);
for k = 1:totalregions
    for i = 2:row-1
        for j = 2:col-1
            val = regionlabel(i,j);
            if val == k
                %%getting all neighbors
                gfneighbor = [gftran(i-1,j-1,k),gftran(i-1,j,k),gftran(i-1,j+1,k),gftran(i,j-1,k)];
                gftran(i,j,k) = min(gfneighbor)+1;
            end
        end
    end
end
%%now doing the second pass
%%decrmenet to go from right top left, bottom to top
for k = 1:totalregions
    for i = row-1:-1:2
        for j = col-1:-1:2
            val = regionlabel(i,j);
            if val == k
                %%getting all neighbors
                gfneighbor = [gftran(i+1,j-1,k),gftran(i+1,j,k),gftran(i+1,j+1,k),gftran(i,j+1,k)];
                gftran(i,j,k) = min(gftran(i,j,k),min(gfneighbor)+1);
            end
        end
    end
end
%%printing grassfire visualization
mysum = sum(gftran,3);
size(mysum);
figure('Name','grassfire transform')
imshow(mysum/max(mysum(:)))

%%