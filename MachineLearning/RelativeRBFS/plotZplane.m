function [] = plotZplane(meanxy,sqsize,Zheight)
%function to plot a a plane in the xy with a height in z

%Specifies vertices
VX = [meanxy(1)+sqsize(1)/2 meanxy(1)-sqsize(1)/2 meanxy(1)-sqsize(1)/2 meanxy(1)+sqsize(1)/2];
VY = [meanxy(2)+sqsize(2)/2 meanxy(2)+sqsize(2)/2 meanxy(2)-sqsize(2)/2 meanxy(2)-sqsize(2)/2];
VZ = [Zheight Zheight Zheight Zheight];

%use patch to plot plane
patch(VX,VY,VZ,'g','FaceAlpha',0.5);


end