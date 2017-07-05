% GIVEN THE COORDINATES (X1, Y1) AND (X2, Y2) OF TWO POINTS

% THIS FUNCTION COMPUTES THE PARAMETERS A, B, C OF THE

% LINE A X + B Y = C THAT IS DEFINED BY THESE TWO POINTS



function [A,B,C]=line_equation(x1,y1,x2,y2)



dif = x1-x2;

% checks for non vertical line segmets

if abs(dif)>0, % non vertical line segment case

   a=(y1-y2)/(x1-x2);

   b=(x1*y2-x2*y1)/(x1-x2);

   A=a;

   B=-1;

   C=-b;

else % vertical line segment case

   A=1;

   B=0;

   C=x1;

end



   