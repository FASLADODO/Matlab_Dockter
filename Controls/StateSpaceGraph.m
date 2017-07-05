function [S] = StateSpaceGraph(b)
m1=10;
m2=350;
kw=500000;
ks=10000;
t = 0:0.01:2;

F=[0,1,0,0;-((ks/m1)+(kw/m1)),-b/m1,ks/m1,b/m1;0,0,1,0;ks/m2,b/m2,-ks/m2,-b/m2];
G=[0;kw/m1;0;0];
H=[1,0,0,0;0,0,1,0];
J=0;

y = step(F,G,H,J,1,t);

plot( t, y(:,1), ':', t, y(:,2), '-' );

S=1;
