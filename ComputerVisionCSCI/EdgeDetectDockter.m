function [thinnedpic,edgepic] = EdgeDetectDockter(picnamestring,definethresh,expansionnum)
%%Rod Dockter
%%Arguments: 
%%picnamestring: picture name string with .jpg
%%edge threshold for detection
%%number of expansions
%%0.38 is good threshold for bw images

%%This function is the edge detection program I develeoped for Programming
%%assignment # 1 which I will reuse for the Hough Transform.
%%This particular edge detecter will make use of multiple operators so as
%%to read in diagonal and straight lines.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Preprocessing to grayscale and image input
%%using double 0 -1 scale
inputpic = imread(picnamestring);
picg = rgb2gray(inputpic);
picg=im2double(picg);
%%picg = Grayscale version
[row,col]=size(picg);
edgepic=zeros(row,col);
edgepic=mat2gray(edgepic);
edgepic=im2double(edgepic);
magthreshold=definethresh;

%%making an array of operator values
opermat = zeros(row,col);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Implementation of sobel and roberts operators
%%Direction and magnitude thresholding occurs in the loop below 

for i=2:row-1
    for j=2:col-1
        %%Produce gradient terms with sobel and roberts operator
        deltas1=picg(i-1,j+1)+2*picg(i,j+1)+picg(i+1,j+1)-picg(i-1,j-1)-2*picg(i,j-1)-picg(i+1,j-1);
        deltas2=picg(i-1,j-1)+2*picg(i-1,j)+picg(i-1,j+1)-picg(i+1,j-1)-2*picg(i+1,j)-picg(i+1,j+1);
        mags=sqrt((deltas1^2)+(deltas2^2));
        opermat(i,j) = mags;
        if mags>magthreshold
            edgepic(i,j)=1;
        end
    end
end
edgepic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Performing expansion on the edge detection image to fill in gaps
%%I will do this four times to be safe with the dual lines produced by the
%%line image.
expansionedge=double(edgepic);
expansionloop=double(edgepic);
for h=1:expansionnum
    for q=2:row-1
       for w=2:col-1
            neighborvector=[expansionloop(q-1,w-1),expansionloop(q-1,w),expansionloop(q-1,w+1),expansionloop(q,w-1),expansionloop(q,w+1),expansionloop(q+1,w-1),expansionloop(q+1,w),expansionloop(q+1,w+1)];
            neighborvectorcount  = 0;
            for k=1:numel(neighborvector)
                if neighborvector(k) == 1
                   neighborvectorcount =1; 
                end
            end
            if neighborvectorcount == 1
                expansionedge(q,w)=1;
            end
       end
    end
    expansionloop=expansionedge;
end

expansionedge=expansionloop;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Edge Thinning of the expanded image using eight connectivity
%%Tracker for changes
change=1;
loopticker = 0;
%%pre allocate
thinnedpic=double(expansionedge);
%% thinnedpicF acts as the un-processed image from which each pass can make decisions about removal
thinnedpicF=double(expansionedge);
while change ~=0
    change=0;
%%%%%%North Pass d:row f:col
    d=row-1;
    f=2;
    while f<col
        pixelN=thinnedpicF(d,f);
        if pixelN == 1
            neighborN=[thinnedpicF(d-1,f-1),thinnedpicF(d-1,f),thinnedpicF(d-1,f+1),thinnedpicF(d,f-1),thinnedpicF(d,f+1),thinnedpicF(d+1,f-1),thinnedpicF(d+1,f),thinnedpicF(d+1,f+1)];
            neighborNcount=0;
            for vc=1:numel(neighborN)
                if neighborN(vc) == 1
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
            elseif neighborN(2) == 1
                trackifN = 1;
            elseif neighborNcount == 2
                trackifN = 1;
            elseif neighborNcount > 2 
                if (neighborN(4)==1) && (neighborN(5)==1)
                    if (neighborN(2)==0) && (neighborN(7)==0)
                        trackifN = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif (neighborN(2)==1) && (neighborN(7)==1)
                    if (neighborN(4)==0) && (neighborN(5)==0)
                        trackifN = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif neighborN(1)==1
                    if (neighborN(2)==0) && (neighborN(4)==0)
                        trackifN = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif neighborN(3)==1
                    if (neighborN(2)==0) && (neighborN(5)==0)
                        trackifN = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif neighborN(6)==1
                    if (neighborN(4)==0) && (neighborN(7)==0)
                        trackifN = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif neighborN(8)==1
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
        if pixelE == 1
            neighborE=[thinnedpicF(d-1,f-1),thinnedpicF(d-1,f),thinnedpicF(d-1,f+1),thinnedpicF(d,f-1),thinnedpicF(d,f+1),thinnedpicF(d+1,f-1),thinnedpicF(d+1,f),thinnedpicF(d+1,f+1)];
            neighborEcount=0;
            for vc=1:numel(neighborE)
                if neighborE(vc) == 1
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
    thinnedpicF=thinnedpic;        elseif neighborE(4) ==1
                trackifE = 1;
            elseif neighborEcount == 2
                trackifE = 1;
            elseif neighborEcount > 2 
                if (neighborE(4)==1) && (neighborE(5)==1)
                    if (neighborE(2)==0) && (neighborE(7)==0)
                        trackifE = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif (neighborE(2)==1) && (neighborE(7)==1)
                    if (neighborE(4)==0) && (neighborE(5)==0)
                        trackifE = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif neighborE(1)==1
                    if (neighborE(2)==0) && (neighborE(4)==0)
                        trackifE = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif neighborE(3)==1
                    if (neighborE(2)==0) && (neighborE(5)==0)
                        trackifE = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif neighborE(6)==1
                    if (neighborE(4)==0) && (neighborE(7)==0)
                        trackifE = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif neighborE(8)==1
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
    thinnedpicF=thinnedpic;
%%%%%%South Pass d:row f:col
    d=2;
    f=2;
    while f<col
        pixelS=thinnedpicF(d,f);
        if pixelS == 1
            neighborS=[thinnedpicF(d-1,f-1),thinnedpicF(d-1,f),thinnedpicF(d-1,f+1),thinnedpicF(d,f-1),thinnedpicF(d,f+1),thinnedpicF(d+1,f-1),thinnedpicF(d+1,f),thinnedpicF(d+1,f+1)];
            neighborScount=0;
            for vc=1:numel(neighborS)
                if neighborS(vc) == 1
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
            elseif neighborS(7) ==1
                trackifS = 1;
            elseif neighborScount == 2
                trackifS = 1;
            elseif neighborScount > 2 
                if (neighborS(4)==1) && (neighborS(5)==1)
                    if (neighborS(2)==0) && (neighborS(7)==0)
                        trackifS = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif (neighborS(2)==1) && (neighborS(7)==1)
                    if (neighborS(4)==0) && (neighborS(5)==0)
                        trackifS = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif neighborS(1)==1
                    if (neighborS(2)==0) && (neighborS(4)==0)
                        trackifS = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif neighborS(3)==1
                    if (neighborS(2)==0) && (neighborS(5)==0)
                        trackifS = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif neighborS(6)==1
                    if (neighborS(4)==0) && (neighborS(7)==0)
                        trackifS = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif neighborS(8)==1
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
    thinnedpicF=thinnedpic;
%%%%%%West Pass d:row f:col
    d=2;
    f=col-1;
    while d<row
        pixelW=thinnedpic(d,f);
        if pixelW == 1
            neighborW=[thinnedpicF(d-1,f-1),thinnedpicF(d-1,f),thinnedpicF(d-1,f+1),thinnedpicF(d,f-1),thinnedpicF(d,f+1),thinnedpicF(d+1,f-1),thinnedpicF(d+1,f),thinnedpicF(d+1,f+1)];
            neighborWcount=0;
            for vc=1:numel(neighborW)
                if neighborW(vc) == 1
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
            elseif neighborW(5) ==1
                trackifW = 1;
            elseif neighborWcount == 2
                trackifW = 1;
            elseif neighborWcount > 2 
                if (neighborW(4)==1) && (neighborW(5)==1)
                    if (neighborW(2)==0) && (neighborW(7)==0)
                        trackifW = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif (neighborW(2)==1) && (neighborW(7)==1)
                    if (neighborW(4)==0) && (neighborW(5)==0)
                        trackifW = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif neighborW(1)==1
                    if (neighborW(2)==0) && (neighborW(4)==0)
                        trackifW = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif neighborW(3)==1
                    if (neighborW(2)==0) && (neighborW(5)==0)
                        trackifW = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif neighborW(6)==1
                    if (neighborW(4)==0) && (neighborW(7)==0)
                        trackifW = 1;
                    else
                        thinnedpic(d,f)=0;
                        change =1;
                    end
                elseif neighborW(8)==1
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