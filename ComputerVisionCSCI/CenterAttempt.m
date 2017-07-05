%%region labeling ellipse
%%
imagestring='rodellipse.jpg';
xcbin = 1;
ycbin = 1;
alphabin = 1;
betabin = 1;
%%This function will work by first attempting to find all the possible
%%centers of ellipses, using a voting scheme on those, then performing a
%%second accumulator array to find the alpha and beta axis lengths
%%I will use the standard form
%%(x-xc)^2/alpha^2+(y-yc)^2/beta^2=1

%%Calling edge detection and thinning algorithm
%%from previous tests I know 0.38 to be a pretty good threshold
[thinnedpic,edgepic] = EdgeDetectDockter(imagestring,0.38,1);
%%
figure('Name','edge image');
imshow(thinnedpic)
%%

%%getting size of image
[row,col]=size(thinnedpic);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%finding center

%%creating arrays of edge points
%%creating arrays of edge points
xcoord = [];
ycoord = [];
index = 1;
for i=1:row
    for j=1:col
        if thinnedpic(i,j) == 1
            xcoord(index) = j;
            ycoord(index) = row - i;
            index=index+1;
        end
    end
end

%% 

length(xcoord)
length(ycoord)

%%

centeracc = zeros(row,col);

randind1 = randi(length(xcoord),1,100*length(xcoord));
randind2 = randi(length(xcoord),1,100*length(xcoord));
randind3 = randi(length(xcoord),1,100*length(xcoord));

for i = 1:100*length(xcoord)
    X1 = xcoord(randind1(i));
    Y1 = ycoord(randind1(i));
    X2 = xcoord(randind2(i));
    Y2 = ycoord(randind2(i));
    X3 = xcoord(randind3(i));
    Y3 = ycoord(randind3(i));
    dist1=sqrt((X1-X2)^2+(Y1-Y2)^2);
    dist2=sqrt((X2-X3)^2+(Y2-Y3)^2);
    dist3=sqrt((X1-X3)^2+(Y1-Y3)^2);
    if (dist1 > 3) && (dist2 > 3) && (dist3 > 3) && (X1 ~= X2) && (X1 ~= X3)
        neighbori = -3:3;
        neighborj = -3:3;
        ntick1=1;
        ntick2=1;
        ntick3=1;
        for ki = neighbori
            for kj = neighborj
                tempi1 = X1+ki;
                tempj1 = Y1+kj;
                tempi2 = X2+ki;
                tempj2 = Y2+kj;
                tempi3 = X3+ki;
                tempj3 = Y3+kj;
                if( thinnedpic(row-tempj1,tempi1) == 1)
                    Xpoints1(ntick1) = tempi1;
                    Ypoints1(ntick1) = tempj1;
                    ntick1=ntick1+1;
                end
                if( thinnedpic(row-tempj2,tempi2) == 1)
                    Xpoints2(ntick2) = tempi2;
                    Ypoints2(ntick2) = tempj2;
                    ntick2=ntick2+1;
                end
                if( thinnedpic(row-tempj3,tempi3) == 1)
                    Xpoints3(ntick3) = tempi3;
                    Ypoints3(ntick3) = tempj3;
                    ntick3=ntick3+1;
                end
            end
        end
    end
    [slope1, intercept1] = LeastSquaresFit(Xpoints1, Ypoints1);
    [slope2, intercept2] = LeastSquaresFit(Xpoints2, Ypoints2);
    [slope3, intercept3] = LeastSquaresFit(Xpoints3, Ypoints3);
    %%calculating intersections
    xint12 = (intercept2 - intercept1)/(slope1-slope2);
    xint13 = (intercept3 - intercept1)/(slope1-slope3);
    yint12 = slope2*xint12 + intercept2;
    yint13 = slope3*xint13 + intercept3;
    %%midpoints between ellipse edge points
    mdp12x = (X1+X2)/2; mdp12y = (Y1+Y2)/2;
    mdp13x = (X1+X3)/2; mdp13y = (Y1+Y3)/2;
    [tslope12, tintercept12] = LeastSquaresFit([xint12,mdp12x], [yint12,mdp12y]);
    [tslope13, tintercept13] = LeastSquaresFit([xint13,mdp13x], [yint13,mdp13y]);
    xcp = (tintercept13 - tintercept12)/(tslope13-tslope12);
    ycp = tslope12*xcp + tintercept12;
    if (round(xcp) < col) && (round(ycp) < row) && (round(xcp) > 0) && (round(ycp) > 0)
       centeracc(round(ycp),round(xcp))=centeracc(round(ycp),round(xcp))+1;
    end
end
%%

[max_vote, mvind] = max(centeracc(:));
[rowi,colj] = ind2sub(size(centeracc),mvind);
rowi
colj
row-rowi

readi = imread('rodellipse.jpg');
readi=im2double(readi);
imshow(readi)
hold on
plot(row-rowi,colj)
hold off








%%
