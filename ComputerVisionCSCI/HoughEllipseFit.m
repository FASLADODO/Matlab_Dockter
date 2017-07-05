%%function [xcparam, ycparam, alphaparam, betaparam] = HoughEllipseFit(imagestring, votethreshold, xcbin, ycbin, alphabin, betabin)
%%Rod Dockter
%%Programming Assignment 2, Part 2
%%Hough transform ellipse fitting

%%The program is written in cell mode because It is easier to run in
%%logical chunks than all at once.
%%If you prefer to run it as its own function just uncomment the function
%%line at the top and end at the bottom
%%
%%Reading image and doing edge detection
imagestring='ellipse1.jpg';
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
imshow(thinnedpic)
%%
%%Finding the size of image and the ranges for edge points
%%getting size of image
[row,col]=size(thinnedpic);

%%creating arrays of edge points
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

%% 

disp('Entering centerfinding')
%%The next section is an algorithm I made up to find the center of the
%%ellipses (roughly)

%%finding the bounds of my edge points so I don't have to check the whole
%%image
xmin = min(xcoord);
xmax = max(xcoord);
ymin=min(ycoord);
ymax=max(ycoord);

%%pre allocating a flagging matrix so I can check edge points of the
%%ellipses
flagmat = zeros(row,col);

%%It works by searching through all the edge points in my image and as soon
%%as it finds an edge point which happens to also have a neighbor to the
%%north or south I begin to examine that line further (hoping it is a line
%%lying on the left side of an ellipse)
%%The algorithm then goes down and up simoultaneously from the edge point
%%untill it no longer sees an edge pixel farther up(or down). Once it
%%reaches the top and bottom it checks to see if the neighbors to left of
%%the top and bottom are zero, the neighbors to NE of the top and the SE of
%%the bottom are zero and that the neighbor to the NW of the top and SW of
%%the bottom is one (an edge) so that I know I am at the left of an edge
%%outline that will continue up and down to the right of me.
%%Once I am assured that the edge I have found is on the left side of an
%%ellipse I calculate the rough ycenter coordinate as the average of the
%%top and bottom of the horizontal line
%%then I move over to the right one step at a time until I find the other
%%side of the ellipse (detect an edge pixel)
%%I then record the distance between the right side and the left side and
%%record the average of the two sides positions as the x center coordinate.
%%I also flag any potential left side as used in a seprate matrix so I
%%don't check the edge again.
%%*Notice that this algorithm is not particularly robust and completely
%%relies on the edges being one pixel thick, the ellipse edges having no
%%gaps, and the ellipse edge image having no pixel inside the ellipse.
%%However I had to find a separate way of determining the centers because
%%my computer did not have enough memory to produce a 4 dimensional
%%accumulator array of sufficient density to search through all xenter and
%%alpha and beta points.

xcenterarray = [];
ycenterarray = [];
cindex = 1;
for i=ymin:ymax
    for j=xmin:xmax
        %%If i am at edge points
        if thinnedpic(i,j) == 1
            downmatch = 0;
            upmatch = 0;
            %%If edge points has not been previously checked
            if flagmat(i,j) == 0
                %%if the point has neighbor above or below it
                if thinnedpic(i-1,j) == 1 || thinnedpic(i+1,j)==1
                    horcoord = j;
                    %%running through edge points above
                    if thinnedpic(i-1,j) == 1
                       newup = i-2;
                       upif = 0;
                       while upif ~= 1
                           %%if the next point up is still an edge
                           %%increment
                           if thinnedpic(newup,horcoord) == 1
                              newup = newup - 1; 
                           elseif thinnedpic(newup,horcoord) == 0
                               %%if top has been found, check to see if
                               %%neighbors satisfy requirements
                               if (thinnedpic(newup,horcoord-1)==0) && (thinnedpic(newup+1,horcoord-1) == 0) && (thinnedpic(newup,horcoord+1) == 1) && (thinnedpic(newup+1,horcoord+1) == 0) 
                                   upmatch = 1;
                                   upmax = newup + 1;
                               end
                               upif = 1;
                           end
                       end
                    end
                    %%running through edge points below
                    if thinnedpic(i+1,j) == 1
                       newdown = i+2;
                       downif = 0;
                       while downif ~= 1
                           %%if next pixel down is an edge, increment
                           if thinnedpic(newdown,horcoord) == 1
                              newdown = newdown + 1; 
                           elseif thinnedpic(newdown,horcoord) == 0
                               %%if bottom has been found, check to see if
                               %%neighbors satisfy requirements
                               if (thinnedpic(newdown,horcoord-1)==0) && (thinnedpic(newdown-1,horcoord-1) == 0) && (thinnedpic(newdown,horcoord+1) == 1) && (thinnedpic(newdown-1,horcoord+1) == 0)
                                   downmatch = 1;
                                   downmax = newdown - 1;
                               end
                               downif = 1;
                           end
                       end
                    end
                end
            end
            %%if the parameters for the top and bottom of edge were met, I
            %%use the horizontal edge as my ellipses left side
            if downmatch == 1 && upmatch == 1
                disp('Found another one')
                fpoints = upmax:downmax;
                %%flagging points
                flagmat(fpoints,horcoord)=1;
                %%Setting ycenter to average of top and bottom
                ycenter = round((upmax+downmax)/2);
                foundside = 0;
                notfail = 0;
                %%Incerementing to the right to find right side of ellipse
                movepnt = horcoord + 20;
                while foundside ~= 1
                    if (thinnedpic(ycenter,movepnt) == 0) && (movepnt <= xmax + 1)
                        movepnt = movepnt + 1;
                    elseif (thinnedpic(ycenter,movepnt) == 1)
                        %%once found I save the column number of the right
                        %%side
                        foundside = 1;
                        rightside = movepnt;
                    else 
                       %%If neither of the above conditions are true I
                       %%break out so I avoid prepetual while loops
                       foundside = 1;
                       notfail = 1;
                    end
                end
                if notfail == 0
                    %%saving the x-center coordinate for each as the
                    %%average of the right and left side column values.
                    xcenter = round((rightside + horcoord)/2);
                    xcenterarray(cindex) = xcenter;
                    ycenterarray(cindex) = ycenter;
                    cindex = cindex + 1;
                end
            end
        end
    end
end

%%Printing values of centers so I can check them
xcenterarray
ycenterarray

%% 
%%Plotting centers to visualize
figure('Name','Ellipse edges with centers');
imshow(thinnedpic)
hold on
plot(xcenterarray,ycenterarray,'.w')
hold off
%%number of ellipses will be..
ellipsenum = length(xcenterarray);

%%

%%now I will make the assumption that that largest and smallest edge pixels
%%will act as the bounds for xcenter and ycenter
%%I will also assume that alpha and beta won't be larger than half the
%%distance between the largest and smallest x and y edge pixels
%%Also beta has to be less than alpha
alphamax = (xmax-xmin)/2;
betamaxtemp = (ymax-ymin)/3;

%%Initiliazing alpha and beta axis values
%%Since I will assume the ellipses all lie fully in the image boundaries
%%then I can limit the range of possible axis value significantly
%%We also have the additional constraint of alpha >= Beta (deal with later)
alphapoints = 1:alphabin:alphamax;
betapoints = 1:betabin:betamax;

%%defining and epsilon since the standard ellipse equation might not equal
%%one always, but something close to it
epsilon = 0.00001;

%%initializing the 3 dimensional accumulator array for centers and axis lengths
%%The first dimension will be the number of center points I found
%%the second and third dimension will be alpha and beta values
Hacc=zeros(ellipsenum,length(alphapoints),length(betapoints));

%%

disp('entering accumulator')
%%Now I will loop through my previously determined center coordinates and
%%test all the possible alpha and beta values for each edge point
%%(x-xc)^2/alpha^2+(y-yc)^2/beta^2=1
%%This allows me to constrain alpha and beta as not being to large
%%With predetermining the center coordinates I can widdle the looping down
%%to 4 for loops
for ellipind = 1:ellipsenum
    ellipind
    for coordind = 1:length(xcoord)
        %%we want explore all possible alpha and beta axis lengths
        %%I will assume all ellipses lie fully in the plain
        dist2c = sqrt((xcoord(coordind)-xcenterarray(ellipind))^2 + (ycoord(coordind)-ycenterarray(ellipind))^2 );
        if ( dist2c <= alphamax)
            aindex = 0;
            for al = alphapoints
                aindex=aindex+1;
                bindex = 0;
                for be = betapoints
                    bindex=bindex+1;
                    %%solving for equation of ellipse
                    ellipse = (((xcoord(coordind)-xcenterarray(ellipind))^2)/(al^2)) + (((ycoord(coordind)-ycenterarray(ellipind))^2)/(be^2));
                    if (abs(ellipse - 1) < epsilon) && (be <=al) && (xcenterarray(ellipind) - al >= xmin-1) && (ycenterarray(ellipind) - be >= ymin-1) %%if equation is satisfied
                        %%incrementing the accumulator array
                        Hacc(ellipind,aindex,bindex) = Hacc(ellipind,aindex,bindex) + 1;
                    end
                end
            end
        end
    end
end

%%
disp('entering vote count')

%%This section will run through the accumulator array and for each center
%%coordinate, pick out the highest voted alpha and beta parameters
%%preallocate
alphaparam = [];
betaparam = [];
%%For each center point (each ellipse) I find the alpha and beta parameters
%%with the highest votes and keep those values
filltick = 1;
for i=1:ellipsenum
    d2mat = squeeze(Hacc(i,:,:));
    [m,ix]=max(d2mat(:));
    [a,b] = ind2sub(size(d2mat),ix);
    alphaparam(filltick)=alphapoints(a);
    betaparam(filltick)=betapoints(b);
    filltick=filltick+1;
end
%% 
%%Code for printing out accumulator arrays (1 at a time)
figure('Name','mesh');
B = squeeze(Hacc(4,:,:));
size(B);
mesh(B)

alphaparam
betaparam
%%

%%Plotting the detected ellipses over the image
figure('Name','Image with Hough Lines');
readin = imread('ellipse1.jpg');
picg = rgb2gray(readin);
picg=im2double(picg);
imshow(picg)
hold on
for i = 1:ellipsenum
    fun = @(x,y) (((x-xcenterarray(i))^2)/(alphaparam(i)^2)) + (((y-ycenterarray(i))^2)/(betaparam(i)^2)) - 1;
    h = ezplot(fun, [1 800 1 600])
    set(h,'Color','r')
    hold on
end
hold off


%%end