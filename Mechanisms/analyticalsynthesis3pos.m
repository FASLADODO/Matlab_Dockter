%Given
alpha2=-45;
alpha3=9.3;
delta2=-31.19;
delta3=-16.34;
p21=2.798;
p31=3.919;
%Assumed
beta2=342.3;
beta3=324.8;
gamma2=30.9;
gamma3=80.6;

%Matrix elements
A1=cosd(beta2)-1;
B1=sind(beta2);
C1=cosd(alpha2)-1;
D1=sind(alpha2);
E1=p21*cosd(delta2);
F1=p21*sind(delta2);

G1=cosd(beta3)-1;
H1=sind(beta3);
I1=cosd(alpha3)-1;
J1=sind(alpha3);
K1=p31*cosd(delta3);
L1=p31*sind(delta3);

A2=cosd(gamma2)-1;
B2=sind(gamma2);
C2=cosd(alpha2)-1;
D2=sind(alpha2);
E2=p21*cosd(delta2);
F2=p21*sind(delta2);

G2=cosd(gamma3)-1;
H2=sind(gamma3);
I2=cosd(alpha3)-1;
J2=sind(alpha3);
K2=p31*cosd(delta3);
L2=p31*sind(delta3);

%%

matleft = [A1,-B1,C1,-D1;B1,A1,D1,C1;G1,-H1,I1,-J1;H1,G1,J1,I1];
matright = [A2,-B2,C2,-D2;B2,A2,D2,C2;G2,-H2,I2,-J2;H2,G2,J2,I2];

ansleft = [E1;F1;K1;L1];
ansright = [E2;F2;K2;L2];

wz = inv(matleft)*ansleft
us = inv(matright)*ansright

w = sqrt(wz(1)^2+wz(2)^2)
theta = atan2(wz(2),wz(1))

z = sqrt(wz(3)^2+wz(4)^2)
phi = atan2(wz(4),wz(3))

u = sqrt(us(1)^2+us(2)^2)
sigma = atan2(us(2),us(1))

s = sqrt(us(3)^2+us(4)^2)
psi = atan2(us(4),us(3))


