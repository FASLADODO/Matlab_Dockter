format short

theta = 45;
d = 0;
a = 40;
alpha = 0;


Amat = [cosd(theta),-sind(theta)*cosd(alpha),sind(theta)*sind(alpha),a*cosd(theta);
    sind(theta),cosd(theta)*cosd(alpha),-cosd(theta)*sind(alpha),a*sind(theta);
    0,sind(alpha),cosd(alpha),d;
    0,0,0,1]

%%
format short

theta = 90;
d = 100;
a = 0;
alpha = -90;


Amat1 = [cosd(theta),-sind(theta)*cosd(alpha),sind(theta)*sind(alpha),a*cosd(theta);
    sind(theta),cosd(theta)*cosd(alpha),-cosd(theta)*sind(alpha),a*sind(theta);
    0,sind(alpha),cosd(alpha),d;
    0,0,0,1]

theta = -90;
d = 0;
a = 0;
alpha = -90;


Amat2 = [cosd(theta),-sind(theta)*cosd(alpha),sind(theta)*sind(alpha),a*cosd(theta);
    sind(theta),cosd(theta)*cosd(alpha),-cosd(theta)*sind(alpha),a*sind(theta);
    0,sind(alpha),cosd(alpha),d;
    0,0,0,1]

theta = -45;
d = 0;
a = 424.3;
alpha = 0;


Amat3 = [cosd(theta),-sind(theta)*cosd(alpha),sind(theta)*sind(alpha),a*cosd(theta);
    sind(theta),cosd(theta)*cosd(alpha),-cosd(theta)*sind(alpha),a*sind(theta);
    0,sind(alpha),cosd(alpha),d;
    0,0,0,1]

theta = 74.74;
d = 0;
a = 403.1;
alpha = 0;


Amat4 = [cosd(theta),-sind(theta)*cosd(alpha),sind(theta)*sind(alpha),a*cosd(theta);
    sind(theta),cosd(theta)*cosd(alpha),-cosd(theta)*sind(alpha),a*sind(theta);
    0,sind(alpha),cosd(alpha),d;
    0,0,0,1]

T = Amat1*Amat2*Amat3*Amat4

%%

theta = -45;

TG = [cosd(theta),-sind(theta),0,778190;
    sind(theta),cosd(theta),0,480350;
    0,0,1,0;
    0,0,0,1];

TC = [1,0,0,0;
    0,1,0,1.5;
    0,0,1,1;
    0,0,0,1];

TO = [1,0,0,5;
    0,1,0,0;
    0,0,1,-2.5;
    0,0,0,1];

format long g

T = TG*TC*TO