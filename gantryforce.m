%https://en.wikipedia.org/wiki/Leadscrew

d = 70 / 1000; %lead screw diameter (mm)

L=deg2rad(29) %lead angle (assuming acme)

T = 1.2 % torqueN-m stepper motor estimate

F = (T*2) / (d*tan(L))
