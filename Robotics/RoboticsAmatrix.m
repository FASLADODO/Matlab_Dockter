% A Matrices
%Enter elements from the DH table

theta = 0;
d = 0.1;
a = 0;
alpha = -90;


MatA = [cosd(theta),-sind(theta)*cosd(alpha),sind(theta)*sind(alpha),a*cosd(theta);
    sind(theta),cosd(theta)*cosd(alpha),-cosd(theta)*sind(alpha),a*sind(theta);
    0,sind(alpha),cosd(alpha),d;
    0,0,0,1];

MatA

%%

theta = 0;
d = 0.1;
a = 0;
alpha = -90;

Mat1 = [cosd(theta),-sind(theta)*cosd(alpha),sind(theta)*sind(alpha),a*cosd(theta);
    sind(theta),cosd(theta)*cosd(alpha),-cosd(theta)*sind(alpha),a*sind(theta);
    0,sind(alpha),cosd(alpha),d;
    0,0,0,1];

theta = -90;
d = 0;
a = 0;
alpha = -90;

Mat2 = [cosd(theta),-sind(theta)*cosd(alpha),sind(theta)*sind(alpha),a*cosd(theta);
    sind(theta),cosd(theta)*cosd(alpha),-cosd(theta)*sind(alpha),a*sind(theta);
    0,sind(alpha),cosd(alpha),d;
    0,0,0,1];

theta =-45;
d = 0;
a = 0.424;
alpha = 0;

Mat3 = [cosd(theta),-sind(theta)*cosd(alpha),sind(theta)*sind(alpha),a*cosd(theta);
    sind(theta),cosd(theta)*cosd(alpha),-cosd(theta)*sind(alpha),a*sind(theta);
    0,sind(alpha),cosd(alpha),d;
    0,0,0,1];

theta = -74.7;
d = 0.4;
alpha = 0;

Mat4 = [cosd(theta),-sind(theta)*cosd(alpha),sind(theta)*sind(alpha),a*cosd(theta);
    sind(theta),cosd(theta)*cosd(alpha),-cosd(theta)*sind(alpha),a*sind(theta);
    0,sind(alpha),cosd(alpha),d;
    0,0,0,1];

Mat1
Mat2
Mat3
Mat4

Mat1*Mat2*Mat3*Mat4

%%

T1 = [0.7,0.7,0,778190;
    -0.7,0.7,0,480352;
    0,0,1,0;
    0,0,0,1];

T2 = [1,0,0,0;
    0,1,0,1.5;
    0,0,1,1;
    0,0,0,1];

T3 = [1,0,0,5;
    0,1,0,0;
    0,0,1,-2.5;
    0,0,0,1];

format long
T1*T2*T3















