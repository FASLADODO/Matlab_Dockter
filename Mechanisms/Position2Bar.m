%Enter these lengths (a ground, b input ...
a = 4;
b = 2;
c = 3;
d = 5;
t2 = 50; %theta2


%compute
Ax = a*cosd(t2);
Ay = a*sind(t2);
P = ((Ay^2)/(Ax-d)^2) + 1;
S = (a^2 - b^2 + c^2 - d^2)/(2*(Ax-d));
R = (d-S)^2 - c^2;
Q = (2*Ay*(d-S))/(Ax-d);

By = [(-Q+sqrt(Q^2-4*P*R))/(2*P), (-Q-sqrt(Q^2-4*P*R))/(2*P)]

Bx = [S - ((2*Ay*By(1))/(2*(Ax-d))), S - ((2*Ay*By(2))/(2*(Ax-d))) ]


%% choose which Bx and By index
index = 1;

B_x = Bx(index);
B_y = By(index);

t3 = atan2((B_y - Ay),(B_x - Ax)) *(180/pi)

t4 = atan2((B_y),(B_x - d)) *(180/pi)