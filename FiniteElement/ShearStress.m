%Set these
Vmax = 1600;
base = 0.5;
height = 2;
dist2top = 0.25;



I = (base*height^3)/12;
Ac = base*dist2top;
Yc = ((height/2.0) - dist2top) + dist2top/2.0;
Q = Ac*Yc;

ShearStress = (Vmax*Q)/(I*base)

%%

s = [2400,2250,1800,1050];
d = [0,0.25,0.5,0.75];

plot(d,s)