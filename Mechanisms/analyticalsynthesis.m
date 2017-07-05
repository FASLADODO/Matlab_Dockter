A = -0.923;
B = 0.997;
C=-0.269;
D = 0.682;
E = -2.318;
F = 0.621;

Zx = -0.249;
Zy = 0.999;


wx = (A*(-C*Zx+D*Zy+E)+B*(-C*Zy-D*Zx+F))/(-2*A)

wy = (A*(-C*Zy-D*Zx+F)+B*(C*Zx-D*Zy-E))/(-2*A)

w = sqrt(wx^2+wy^2)
theta = atand(wy/wx)