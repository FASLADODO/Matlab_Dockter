function [maggrad,direcgrad,edgepic,thinnedpic,expansionedge] = EdgeDetectadjust(picnamestring,operatortype,gradientdirection,definethresh)
%%Arguments: 
%%picnamestring: picture name string with .jpg
%%operatortype: 'roberts' or 'sobel'
%%directionthreshold 'x' or 'y' or 'xy' (for both)
if nargin == 1
   operatortype = 'sobel';
   gradientdirection = 'xy';
end
if nargin == 2
   gradientdirection='xy';
end

%%Rod Dockter Feb 2012
%%CSCI 5561 Prof. Papanikolopolous
%%Programming Assignment 1
%%URL: https://sites.google.com/a/umn.edu/umnworkdockter/edgedetection
%%This programmed is designed to read in an image, convert it to grayscale
%%Then perform statistical threshold calculations on the image.
%%Then using either the Roberts or Sobel(3x3) operator, detect edges 
%%Which correspond to pixel gradients which exceed this threshold.
%%After that, the program will perform expansion to fill gaps.
%%Then finally Thin the connectivity to remove overly thick edges.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Preprocessing to grayscale and image input
inputpic = imread(picnamestring);
picg = rgb2gray(inputpic);
picg=im2uint8(picg);
save picg
%%picg = Grayscale version
[row,col]=size(picg);
%%magnitude gradient
maggrad=zeros(row,col);
maggrad=mat2gray(maggrad);
maggrad=im2uint8(maggrad);
save maggrad
%% direction gradient
direcgrad=zeros(row,col);
direcgrad=mat2gray(direcgrad);
direcgrad=im2uint8(direcgrad);
save direcgrad
edgepic=zeros(row,col);
edgepic=mat2gray(edgepic);
edgepic=im2uint8(edgepic);
save edgepic

%%Defining relevant operators
%Roberts
%%R1=[0,1;-1,0];
%%R2=[1,0;0,-1];
%%Sobel
%%S1=[-1,0,1;-2,0,2;-1,0,1];
%%S2=[1,2,1;0,0,0;-1,-2,-1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Determining Threshold
%%Using Otsu Method of between class variance for magnitude threshold
%%maxgray=max(picg(:));
totalpixel=row*col;
%%turns the picture into a vector to perform a histogram on the whole
%%matrix (1:256 since arrays are based at one
grayhist=zeros(1,256);
picvector=reshape(picg,1,totalpixel);
for val=1:totalpixel
    histvalue=round(picvector(val))+1;
    grayhist(histvalue)=grayhist(histvalue)+1;
end 
%%Perform loop to determine weights and means for the forground and
%%background of the histogram
%%*I consulted wikipedia for specifics of this method

%%Initialize values
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
    weightedtotal=(k-1)*grayhist(k);
end
%%Loop for threshold determination (Maximizing between class variance)
for m=1:1:256
    Wbnum=Wbnum+grayhist(m);
    Wfnum=totalpixel-Wbnum;
    Wb=Wbnum/Wbden;
    Wf=Wfnum/Wfden;
    if (Wfnum ~= 0) && (Wbnum~=0)
        Ubnum=Ubnum + m*grayhist(m);
        Ubden=Wbnum;
        Ub=Ubnum/Ubden;
        Ufnum=weightedtotal-Ubnum;
        Ufden=Wfnum;
        Uf=Ufnum/Ufden;
    end
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
if nargin == 4
magthreshold=definethresh;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Implementation of sobel or roberts operators
%%Direction and magnitude thresholding occurs in the loop below 
%%Case structure is used for direction threshold
switch operatortype
    case 'roberts'
        %%row represents y from top to bottom
        %%col represents x from left to right
        for i=1:row-1
            for j=1:col-1
            %%produce gradient terms with roberts operator
            deltar1=picg(i,j+1)-picg(i+1,j);
            deltar2=picg(i,j)-picg(i+1,j+1);
            magr=uint8(sqrt(double(deltar1)^2+double(deltar2)^2));
            directionr=uint8(atan2(double(deltar2),double(deltar1)));
            maggrad(i,j)=round(magr);
            direcgrad(i,j)=uint8((directionr+pi)*(255/(2*pi)));
            switch gradientdirection
                case 'x'
                    if magr>magthreshold
                        if (abs(directionr)>=(3*pi)/4) || (abs(directionr)<=(pi)/4)
                            edgepic(i,j)=255;
                        end
                    end
                case 'y'
                    if magr>magthreshold
                        if (abs(directionr)<=(3*pi)/4) && (abs(directionr)>=(pi)/4)
                            edgepic(i,j)=255;
                        end
                    end
                case 'xy'
                    if magr>magthreshold
                        edgepic(i,j)=255;
                    end
            end
            end
        end
    case 'sobel'
        for i=2:row-1
            for j=2:col-1
            %%Produce gradient terms with sobel operator
            deltas1=picg(i-1,j+1)+uint8(2)*picg(i,j+1)+picg(i+1,j+1)-picg(i-1,j-1)-uint8(2)*picg(i,j-1)-picg(i+1,j-1);
            deltas2=picg(i-1,j-1)+uint8(2)*picg(i-1,j)+picg(i-1,j+1)-picg(i+1,j-1)-uint8(2)*picg(i+1,j)-picg(i+1,j+1);
            mags=uint8(sqrt(double(deltas1)^2+double(deltas2)^2));
            directions=uint8(atan2(double(deltas2),double(deltas1)));
            maggrad(i,j)=round(mags);
            direcgrad(i,j)=uint8((directions+pi)*(255/(2*pi)));
            switch gradientdirection
                case 'x'
                    if mags>magthreshold
                        if (abs(directions)>=(3*pi)/4) || (abs(directions)<=(pi)/4)
                            edgepic(i,j)=255;
                        end
                    end
                case 'y'
                    if mags>magthreshold
                        if (directions<=(3*pi)/4) && (directions>=(pi)/4)
                            edgepic(i,j)=255;
                        end
                    end
                case 'xy'
                    if mags>magthreshold
                        edgepic(i,j)=255;
                    end
            end
            end
        end
end
maggrad;
direcgrad;
edgepic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Performing expansion on the edge detection image to fill in gaps
expansionedge=uint8(edgepic);
for q=2:row-1
   for w=2:col-1
        neighborvector=[edgepic(q-1,w-1),edgepic(q-1,w),edgepic(q-1,w+1),edgepic(q,w-1),edgepic(q,w+1),edgepic(q+1,w-1),edgepic(q+1,w),edgepic(q+1,w+1)];
        neighborvectorcount  = 0;
        for vc=1:numel(neighborvector)
            if neighborvector(vc) == 255
               neighborvectorcount =1; 
            end
        end
        if neighborvectorcount == 1
            expansionedge(q,w)=255;
        end
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Edge Thinning of the expanded image using eight connectivity
%%Tracker for changes
change=1;
loopticker = 0;
%%pre allocate
thinnedpic=uint8(expansionedge);
%% thinnedpicF acts as the un-processed image from which each pass can make decisions about removal
thinnedpicF=uint8(expansionedge);
while change ~=0
    change=0;
%%%%%%North Pass d:row f:col
    d=row-1;
    f=2;
    while f<col
        pixelN=thinnedpicF(d,f);
        if pixelN == 255
            neighborN=[thinnedpicF(d-1,f-1),thinnedpicF(d-1,f),thinnedpicF(d-1,f+1),thinnedpicF(d,f-1),thinnedpicF(d,f+1),thinnedpicF(d+1,f-1),thinnedpicF(d+1,f),thinnedpicF(d+1,f+1)];
            neighborNcount=0;
            for vc=1:numel(neighborN)
                if neighborN(vc) == 255
                    neighborNcount = neighborNcount + 1;
                end
            end
            trackifN = 0;
            if neighborNcount <= 1
                trackifN = 1;
            elseif neighborNcount == 8
                trackifN = 1;
            elseif neighborNcount == 7
                trackifN = 1;                
            elseif neighborN(2) ==255
                trackifN = 1;
            elseif neighborNcount == 2
                if (neighborN(1)==255) && (neighborN(2)==255)
                    thinnedpic(d,f)=0;
                    change =1;
                elseif (neighborN(2)==255) && (neighborN(3)==255)
                    thinnedpic(d,f)=0;
                    change =1;
                elseif (neighborN(1)==255) && (neighborN(4)==255)
                    thinnedpic(d,f)=0;
                    change =1;
                elseif (neighborN(3)==255) && (neighborN(5)==255)
                    thinnedpic(d,f)=0;
                    change =1;
                elseif (neighborN(4)==255) && (neighborN(6)==255)
                    thinnedpic(d,f)=0;
                    change =1;
                elseif (neighborN(6)==255) && (neighborN(7)==255)
                    thinnedpic(d,f)=0;
                    change =1;
                elseif (neighborN(7)==255) && (neighborN(8)==255)
                    thinnedpic(d,f)=0;
                    change =1;
                elseif (neighborN(8)==255) && (neighborN(5)==255)
                    thinnedpic(d,f)=0;
                    change =1;
                else
                    trackifN=1;
                end
            elseif neighborNcount > 2 
                if (neighborN(4)==255) && (neighborN(5)==255)
                    if (neighborN(2)==0) && (neighborN(7)==0)
                        trackifN = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif (neighborN(2)==255) && (neighborN(7)==255)
                    if (neighborN(4)==0) && (neighborN(5)==0)
                        trackifN = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif neighborN(1)==255
                    if (neighborN(2)==0) && (neighborN(4)==0)
                        trackifN = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif neighborN(3)==255
                    if (neighborN(2)==0) && (neighborN(5)==0)
                        trackifN = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif neighborN(6)==255
                    if (neighborN(4)==0) && (neighborN(7)==0)
                        trackifN = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif neighborN(8)==255
                    if (neighborN(5)==0) && (neighborN(7)==0)
                        trackifN = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                else
                    thinnedpic(d,f)=0;
                    change =1;
                end
            elseif trackifN == 0
                thinnedpic(d,f)=0;
                change=1;
            end
        end
        d=d-1;
        if d==1
            d=row-1;
            f=f+1;
        end
    end
%%%%%%East Pass d:row f:col
    d=2;
    f=2;
    while d<row
        pixelE=thinnedpicF(d,f);
        if pixelE == 255
            neighborE=[thinnedpicF(d-1,f-1),thinnedpicF(d-1,f),thinnedpicF(d-1,f+1),thinnedpicF(d,f-1),thinnedpicF(d,f+1),thinnedpicF(d+1,f-1),thinnedpicF(d+1,f),thinnedpicF(d+1,f+1)];
            neighborEcount=0;
            for vc=1:numel(neighborE)
                if neighborE(vc) == 255
                    neighborEcount = neighborEcount + 1;
                end
            end
            trackifE = 0;
            if neighborEcount <= 1
                trackifE = 1;
            elseif neighborEcount == 8
                trackifE = 1;
            elseif neighborEcount == 7
                trackifE=1;
            elseif neighborE(4) ==255
                trackifE = 1;
            elseif neighborEcount == 2
                if (neighborE(1)==255) && (neighborE(2)==255)
                    thinnedpic(d,f)=0;
                    change =1;
                elseif (neighborE(2)==255) && (neighborE(3)==255)
                    thinnedpic(d,f)=0;
                    change =1;
                elseif (neighborE(1)==255) && (neighborE(4)==255)
                    thinnedpic(d,f)=0;
                    change =1;
                elseif (neighborE(3)==255) && (neighborE(5)==255)
                    thinnedpic(d,f)=0;
                    change =1;
                elseif (neighborE(4)==255) && (neighborE(6)==255)
                    thinnedpic(d,f)=0;
                    change =1;
                elseif (neighborE(6)==255) && (neighborE(7)==255)
                    thinnedpic(d,f)=0;
                    change =1;
                elseif (neighborE(7)==255) && (neighborE(8)==255)
                    thinnedpic(d,f)=0;
                    change =1;
                elseif (neighborE(8)==255) && (neighborE(5)==255)
                    thinnedpic(d,f)=0;
                    change =1;
                else
                    trackifE=1;
                end
            elseif neighborEcount > 2 
                if (neighborE(4)==255) && (neighborE(5)==255)
                    if (neighborE(2)==0) && (neighborE(7)==0)
                        trackifE = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif (neighborE(2)==255) && (neighborE(7)==255)
                    if (neighborE(4)==0) && (neighborE(5)==0)
                        trackifE = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif neighborE(1)==255
                    if (neighborE(2)==0) && (neighborE(4)==0)
                        trackifE = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif neighborE(3)==255
                    if (neighborE(2)==0) && (neighborE(5)==0)
                        trackifE = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif neighborE(6)==255
                    if (neighborE(4)==0) && (neighborE(7)==0)
                        trackifE = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif neighborE(8)==255
                    if (neighborE(5)==0) && (neighborE(7)==0)
                        trackifE = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                else
                    thinnedpic(d,f)=0;
                    change =1;
                end
            elseif trackifE == 0
                thinnedpic(d,f)=0;
                change=1;
            end
        end
        f=f+1;
        if f==col
            f=2;
            d=d+1;
        end
    end
%%%%%%South Pass d:row f:col
    d=2;
    f=2;
    while f<col
        pixelS=thinnedpicF(d,f);
        if pixelS == 255
            neighborS=[thinnedpicF(d-1,f-1),thinnedpicF(d-1,f),thinnedpicF(d-1,f+1),thinnedpicF(d,f-1),thinnedpicF(d,f+1),thinnedpicF(d+1,f-1),thinnedpicF(d+1,f),thinnedpicF(d+1,f+1)];
            neighborScount=0;
            for vc=1:numel(neighborS)
                if neighborS(vc) == 255
                    neighborScount = neighborScount + 1;
                end
            end
            trackifS = 0;
            if neighborScount <= 1
                trackifS = 1;
            elseif neighborScount == 8
                trackifS = 1; 
            elseif neighborScount == 7
                trackifS=1;
            elseif neighborS(7) ==255
                trackifS = 1;
            elseif neighborScount == 2
                if (neighborS(1)==255) && (neighborS(2)==255)
                    thinnedpic(d,f)=0;
                    change =1;
                elseif (neighborS(2)==255) && (neighborS(3)==255)
                    thinnedpic(d,f)=0;
                    change =1;
                elseif (neighborS(1)==255) && (neighborS(4)==255)
                    thinnedpic(d,f)=0;
                    change =1;
                elseif (neighborS(3)==255) && (neighborS(5)==255)
                    thinnedpic(d,f)=0;
                    change =1;
                elseif (neighborS(4)==255) && (neighborS(6)==255)
                    thinnedpic(d,f)=0;
                    change =1;
                elseif (neighborS(6)==255) && (neighborS(7)==255)
                    thinnedpic(d,f)=0;
                    change =1;
                elseif (neighborS(7)==255) && (neighborS(8)==255)
                    thinnedpic(d,f)=0;
                    change =1;
                elseif (neighborS(8)==255) && (neighborS(5)==255)
                    thinnedpic(d,f)=0;
                    change =1;
                else
                    trackifS=1;
                end
            elseif neighborScount > 2 
                if (neighborS(4)==255) && (neighborS(5)==255)
                    if (neighborS(2)==0) && (neighborS(7)==0)
                        trackifS = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif (neighborS(2)==255) && (neighborS(7)==255)
                    if (neighborS(4)==0) && (neighborS(5)==0)
                        trackifS = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif neighborS(1)==255
                    if (neighborS(2)==0) && (neighborS(4)==0)
                        trackifS = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif neighborS(3)==255
                    if (neighborS(2)==0) && (neighborS(5)==0)
                        trackifS = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif neighborS(6)==255
                    if (neighborS(4)==0) && (neighborS(7)==0)
                        trackifS = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif neighborS(8)==255
                    if (neighborS(5)==0) && (neighborS(7)==0)
                        trackifS = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                else
                    thinnedpic(d,f)=0;
                    change =1;
                end
            elseif trackifS == 0
                thinnedpic(d,f)=0;
                change=1;
            end
        end
        d=d+1;
        if d==row
            d=2;
            f=f+1;
        end
    end
    
%%%%%%West Pass d:row f:col
    d=2;
    f=col-1;
    while d<row
        pixelW=thinnedpic(d,f);
        if pixelW == 255
            neighborW=[thinnedpicF(d-1,f-1),thinnedpicF(d-1,f),thinnedpicF(d-1,f+1),thinnedpicF(d,f-1),thinnedpicF(d,f+1),thinnedpicF(d+1,f-1),thinnedpicF(d+1,f),thinnedpicF(d+1,f+1)];
            neighborWcount=0;
            for vc=1:numel(neighborW)
                if neighborW(vc) == 255
                    neighborWcount = neighborWcount + 1;
                end
            end
            trackifW = 0;
            if neighborWcount <= 1
                trackifW = 1;
            elseif neighborWcount == 8
                trackifW = 1;  
            elseif neighborWcount == 7
                trackifW=1;
            elseif neighborW(5) ==255
                trackifW = 1;
            elseif neighborWcount == 2
                if (neighborW(1)==255) && (neighborW(2)==255)
                    thinnedpic(d,f)=0;
                    change =1;
                elseif (neighborW(2)==255) && (neighborW(3)==255)
                    thinnedpic(d,f)=0;
                    change =1;
                elseif (neighborW(1)==255) && (neighborW(4)==255)
                    thinnedpic(d,f)=0;
                    change =1;
                elseif (neighborW(3)==255) && (neighborW(5)==255)
                    thinnedpic(d,f)=0;
                    change =1;
                elseif (neighborW(4)==255) && (neighborW(6)==255)
                    thinnedpic(d,f)=0;
                    change =1;
                elseif (neighborW(6)==255) && (neighborW(7)==255)
                    thinnedpic(d,f)=0;
                    change =1;
                elseif (neighborW(7)==255) && (neighborW(8)==255)
                    thinnedpic(d,f)=0;
                    change =1;
                elseif (neighborW(8)==255) && (neighborW(5)==255)
                    thinnedpic(d,f)=0;
                    change =1;
                else
                    trackifW=1;
                end
            elseif neighborWcount > 2 
                if (neighborW(4)==255) && (neighborW(5)==255)
                    if (neighborW(2)==0) && (neighborW(7)==0)
                        trackifW = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif (neighborW(2)==255) && (neighborW(7)==255)
                    if (neighborW(4)==0) && (neighborW(5)==0)
                        trackifW = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif neighborW(1)==255
                    if (neighborW(2)==0) && (neighborW(4)==0)
                        trackifW = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif neighborW(3)==255
                    if (neighborW(2)==0) && (neighborW(5)==0)
                        trackifW = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif neighborW(6)==255
                    if (neighborW(4)==0) && (neighborW(7)==0)
                        trackifW = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif neighborW(8)==255
                    if (neighborW(5)==0) && (neighborW(7)==0)
                        trackifW = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                else
                    thinnedpic(d,f)=0;
                    change =1;
                end
            elseif trackifW == 0
                thinnedpic(d,f)=0;
                change=1;
            end
        end
        f=f-1;
        if f==1
            f=col-1;
            d=d+1;
        end
    end
    loopticker=loopticker+1;
    %%update un processed to loop proccessed for next pass
    thinnedpicF=thinnedpic;
end
thinnedpic;
end