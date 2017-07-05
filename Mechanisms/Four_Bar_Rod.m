function [x1,y1,x2,y2,x3,y3,x4,y4] = Four_Bar_Rod(a,b,c,d,theta2)

%Using Norton Pg 153

A_x = a*cosd(theta2);
A_y = a*sind(theta2);

S = (a^2-b^2+c^2-d^2)/(2*(A_x - d));

%Finding B_y
P= ((A_y^2)/((A_x-d)^2)) + 1;
Q = (2*A_y*(d-S))/(A_x-d);
R = (d-S)^2 - c^2;

B_y_options = [(-Q +sqrt(Q^2 - 4 * P * R))/ (2*P) , (-Q - sqrt(Q^2 - 4 * P * R))/(2*P)]

B_x_options = [S - ((2*A_y*B_y_options(1))/(2*(A_x-d))), S - ((2*A_y*B_y_options(2))/(2*(A_x-d)))];

if abs(sqrt(B_x_options(1)^2 + B_y_options(1)^2) - c) < abs(sqrt(B_x_options(2)^2 + B_y_options(2)^2) - c)
    B_x = B_y_options(1);
    B_y = B_y_options(1);
else
    B_x = B_y_options(2);
    B_y = B_y_options(2);
end

theta3 = atand((B_y - A_y)/(B_x-A_x));
theta4 = atand((B_y)/(B_x-d));

x1=0;
y1=0;
x2=A_x;
y2=A_y;
x3=B_x;
y3=B_y;
x4=d;
y4=0;
