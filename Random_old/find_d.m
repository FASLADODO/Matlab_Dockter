function [Y] = find_d(diam)
rad=diam/2.0;
yc=162.65;
xc=-3.2;
N=25;
d=79.544;
y1=-sqrt(rad^2-(N-xc)^2)+yc;
y2=-sqrt(rad^2-(N+d-xc)^2)+yc;
y3=sqrt(rad^2-(N+d-xc)^2)+yc;
y4=sqrt(rad^2-(N-xc)^2)+yc;
Y=[y1,y2,y3,y4];

