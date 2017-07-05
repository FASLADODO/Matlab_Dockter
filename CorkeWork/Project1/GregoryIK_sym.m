function [ q_vec ] = GregoryIK_sym( Transform, jointLengths )
% Symbolic solution for gregory IK
% using output from corke toolbox .ikine_sym()

%choose solution version
sol = 2;

%joint lengths
L1 = jointLengths.L1;
L2 = jointLengths.L2;

%unpack transform data
% nx ox ax tx
% ny oy ay ty
% nz oz az tz
nx = Transform(1,1);
ny = Transform(2,1);
nz = Transform(3,1);
ox = Transform(1,2);
oy = Transform(2,2);
oz = Transform(3,2);
ax = Transform(1,3);
ay = Transform(2,3);
az = Transform(3,3);
tx = Transform(1,4);
ty = Transform(2,4);
tz = Transform(3,4);

%Solution 1
if(sol == 1)
    %Joint 1
    q1 = atan2(-ty, -tx);
    S1 = sin(q1);
    C1 = cos(q1);

    %Joint 2
    q2 = atan2(- C1^2*tx^2 + C1^2*ty^2 - 2*S1*C1*tx*ty - L1^2 + L2^2 - ty^2 - tz^2, -(4*C1^6*L1^2*ty^2 - 12*C1^4*L1^2*ty^2 - C1^4*tx^4 + 6*C1^4*tx^2*ty^2 - C1^4*ty^4 + 8*C1^3*L1^2*S1*tx*ty - 4*C1^3*S1*tx^3*ty + 4*C1^3*S1*tx*ty^3 + 2*C1^2*L1^2*tx^2 + 10*C1^2*L1^2*ty^2 + 2*C1^2*L2^2*tx^2 - 2*C1^2*L2^2*ty^2 - 6*C1^2*tx^2*ty^2 - 2*C1^2*tx^2*tz^2 + 2*C1^2*ty^4 + 2*C1^2*ty^2*tz^2 + 8*C1*L1^2*S1^3*tx*ty - 4*C1*L1^2*S1*tx*ty + 4*C1*L2^2*S1*tx*ty - 4*C1*S1*tx*ty^3 - 4*C1*S1*tx*ty*tz^2 - L1^4 + 2*L1^2*L2^2 + 4*L1^2*S1^6*ty^2 - 2*L1^2*ty^2 + 2*L1^2*tz^2 - L2^4 + 2*L2^2*ty^2 + 2*L2^2*tz^2 - ty^4 - 2*ty^2*tz^2 - tz^4)^(1/2)) - angle(-L1*(ty*C1^2*S1*i + tx*C1*i + ty*S1^3*i - tz));
    S2 = sin(q2);
    C2 = cos(q2);

    %Joint 3
    q3 = atan2(L1 + S2*tz - C1*C2*tx - C2*S1*ty, - C2*tz - C1*S2*tx - S1*S2*ty);
    S3 = sin(q3);
    C3 = cos(q3);

    %Joint 4
    q4 = atan2(S1*ax - C1*ay, C2*S3*az + C3*S2*az - C1*C2*C3*ax - C2*C3*S1*ay + C1*S2*S3*ax + S1*S2*S3*ay);
    S4 = sin(q4);
    C4 = cos(q4);

    %Joint 5
    q5 = atan2(S2*S3*nx - C2*C3*nx - C1*C2*S3*nx - C1*C3*S2*nx - C2*S1*S3*ny - C3*S1*S2*ny, C1*S4*ny - S1*S4*nx - C2*C4*S3*nx - C3*C4*S2*nx - C4*S1*S2*S3*ny + C1*C2*C3*C4*nx + C2*C3*C4*S1*ny - C1*C4*S2*S3*nx);
    S5 = sin(q5);
    C5 = cos(q5);

    %Joint 6
    q6 = atan2(1, -(nx^2 + ny^2 + C1^2*nx^2 + C4^2*nx^2 - C1^2*ny^2 - C4^2*ny^2 - 2*C1^2*C4^2*nx^2 - C1^2*C5^2*nx^2 - C2^2*C4^2*nx^2 + C2^2*C5^2*nx^2 - C3^2*C4^2*nx^2 + C3^2*C5^2*nx^2 - C4^2*C5^2*nx^2 + 2*C1^2*C4^2*ny^2 + C1^2*C5^2*ny^2 + C2^2*C4^2*ny^2 - C2^2*C5^2*ny^2 + C3^2*C4^2*ny^2 - C3^2*C5^2*ny^2 + C4^2*C5^2*ny^2 + C1^2*C2^2*C4^2*nx^2 - C1^2*C2^2*C5^2*nx^2 + C1^2*C3^2*C4^2*nx^2 - C1^2*C3^2*C5^2*nx^2 + 2*C2^2*C3^2*C4^2*nx^2 + 2*C1^2*C4^2*C5^2*nx^2 - 2*C2^2*C3^2*C5^2*nx^2 + C2^2*C4^2*C5^2*nx^2 + C3^2*C4^2*C5^2*nx^2 - C1^2*C2^2*C4^2*ny^2 + C1^2*C2^2*C5^2*ny^2 - C1^2*C3^2*C4^2*ny^2 + C1^2*C3^2*C5^2*ny^2 - 2*C2^2*C3^2*C4^2*ny^2 - 2*C1^2*C4^2*C5^2*ny^2 + 2*C2^2*C3^2*C5^2*ny^2 - C2^2*C4^2*C5^2*ny^2 - C3^2*C4^2*C5^2*ny^2 + 2*C1*S1*nx*ny - 2*C1*C4*C5*S5*nx^2 - 4*C1*C4^2*S1*nx*ny - 2*C1*C5^2*S1*nx*ny - 2*C1*C2*C4^2*S2*nx^2 + 2*C1*C2*C5^2*S2*nx^2 - 2*C1*C3*C4^2*S3*nx^2 + 2*C1*C3*C5^2*S3*nx^2 - 2*C1^2*C2^2*C3^2*C4^2*nx^2 + 2*C1^2*C2^2*C3^2*C5^2*nx^2 - C1^2*C2^2*C4^2*C5^2*nx^2 - C1^2*C3^2*C4^2*C5^2*nx^2 - 2*C2^2*C3^2*C4^2*C5^2*nx^2 + 2*C1^2*C2^2*C3^2*C4^2*ny^2 - 2*C1^2*C2^2*C3^2*C5^2*ny^2 + C1^2*C2^2*C4^2*C5^2*ny^2 + C1^2*C3^2*C4^2*C5^2*ny^2 + 2*C2^2*C3^2*C4^2*C5^2*ny^2 + 2*C1*C2^2*C4^2*S1*nx*ny - 2*C1*C2^2*C5^2*S1*nx*ny + 2*C1*C3^2*C4^2*S1*nx*ny - 2*C1*C3^2*C5^2*S1*nx*ny + 4*C1*C4^2*C5^2*S1*nx*ny + 4*C1*C2*C3^2*C4^2*S2*nx^2 - 4*C1*C2*C3^2*C5^2*S2*nx^2 + 4*C1*C2^2*C3*C4^2*S3*nx^2 + 2*C1*C2*C4^2*C5^2*S2*nx^2 - 4*C1*C2^2*C3*C5^2*S3*nx^2 + 2*C1*C3*C4^2*C5^2*S3*nx^2 + 2*C2*C3*C4*S4*nx*ny + 2*C1^2*C2^2*C3^2*C4^2*C5^2*nx^2 - 2*C1^2*C2^2*C3^2*C4^2*C5^2*ny^2 - 2*C4*C5*S1*S5*nx*ny - 2*C4*S2*S3*S4*nx*ny - 2*C2*C4*C5*S2*S5*nx^2 - 2*C3*C4*C5*S3*S5*nx^2 + 2*C2*C4*C5*S2*S5*ny^2 + 2*C3*C4*C5*S3*S5*ny^2 - 2*C2*C4*S1*S3*S4*nx^2 - 2*C3*C4*S1*S2*S4*nx^2 - 2*C2*C4^2*S1*S2*nx*ny + 2*C2*C5^2*S1*S2*nx*ny - 2*C3*C4^2*S1*S3*nx*ny + 2*C3*C5^2*S1*S3*nx*ny + 4*C1*C2^2*C4*C5*S5*nx^2 + 4*C1*C3^2*C4*C5*S5*nx^2 - 2*C2*C3*C4^2*S2*S3*nx^2 + 2*C2*C3*C5^2*S2*S3*nx^2 + 2*C2*C3*C4^2*S2*S3*ny^2 - 2*C2*C3*C5^2*S2*S3*ny^2 + 4*C2^2*C4*C5*S1*S5*nx*ny + 4*C3^2*C4*C5*S1*S5*nx*ny - 2*C5*S1*S2*S3*S4*S5*nx^2 + 4*C1^2*C4*S2*S3*S4*nx*ny + 2*C4*C5^2*S2*S3*S4*nx*ny + 2*C1^2*C2*C4*C5*S2*S5*nx^2 + 4*C2*C3^2*C4*C5*S2*S5*nx^2 + 2*C1^2*C3*C4*C5*S3*S5*nx^2 + 4*C2^2*C3*C4*C5*S3*S5*nx^2 - 2*C1^2*C2*C4*C5*S2*S5*ny^2 - 4*C2*C3^2*C4*C5*S2*S5*ny^2 - 2*C1^2*C3*C4*C5*S3*S5*ny^2 - 4*C2^2*C3*C4*C5*S3*S5*ny^2 + 2*C2*C4*C5^2*S1*S3*S4*nx^2 + 2*C3*C4*C5^2*S1*S2*S4*nx^2 + 4*C2*C3^2*C4^2*S1*S2*nx*ny - 4*C2*C3^2*C5^2*S1*S2*nx*ny + 4*C2^2*C3*C4^2*S1*S3*nx*ny + 2*C2*C4^2*C5^2*S1*S2*nx*ny - 4*C2^2*C3*C5^2*S1*S3*nx*ny + 2*C3*C4^2*C5^2*S1*S3*nx*ny - 8*C1*C2^2*C3^2*C4*C5*S5*nx^2 + 2*C1^2*C2*C3*C4^2*S2*S3*nx^2 - 2*C1^2*C2*C3*C5^2*S2*S3*nx^2 + 2*C2*C3*C4^2*C5^2*S2*S3*nx^2 - 2*C1^2*C2*C3*C4^2*S2*S3*ny^2 + 2*C1^2*C2*C3*C5^2*S2*S3*ny^2 - 2*C2*C3*C4^2*C5^2*S2*S3*ny^2 + 2*C1*C2*C4*S3*S4*nx*ny + 2*C1*C3*C4*S2*S4*nx*ny - 4*C1*C2^2*C3^2*C4^2*S1*nx*ny + 4*C1*C2^2*C3^2*C5^2*S1*nx*ny - 2*C1*C2^2*C4^2*C5^2*S1*nx*ny - 2*C1*C3^2*C4^2*C5^2*S1*nx*ny + 2*C2*C5*S3*S4*S5*nx*ny + 2*C3*C5*S2*S4*S5*nx*ny + 2*C1*C2*C3*C4*S1*S4*nx^2 - 2*C1*C2*C3*C4*S1*S4*ny^2 - 4*C1*C2*C3^2*C4^2*C5^2*S2*nx^2 - 4*C1*C2^2*C3*C4^2*C5^2*S3*nx^2 + 2*C2*C3*C5*S1*S4*S5*nx^2 - 4*C1^2*C2*C3*C4*S4*nx*ny - 2*C2*C3*C4*C5^2*S4*nx*ny - 2*C1*C4*S1*S2*S3*S4*nx^2 + 2*C1*C4*S1*S2*S3*S4*ny^2 - 4*C2*C3^2*C4^2*C5^2*S1*S2*nx*ny - 4*C2^2*C3*C4^2*C5^2*S1*S3*nx*ny + 2*C1*C5*S2*S3*S4*S5*nx*ny - 2*C1^2*C2*C3*C4^2*C5^2*S2*S3*nx^2 + 2*C1^2*C2*C3*C4^2*C5^2*S2*S3*ny^2 + 2*C1*C2*C5*S1*S3*S4*S5*nx^2 + 2*C1*C3*C5*S1*S2*S4*S5*nx^2 - 2*C1*C2*C5*S1*S3*S4*S5*ny^2 - 2*C1*C3*C5*S1*S2*S4*S5*ny^2 - 2*C1*C2*C4*C5^2*S3*S4*nx*ny - 2*C1*C3*C4*C5^2*S2*S4*nx*ny + 4*C1*C2^2*C3^2*C4^2*C5^2*S1*nx*ny - 4*C1^2*C2*C5*S3*S4*S5*nx*ny - 4*C1^2*C3*C5*S2*S4*S5*nx*ny - 2*C1*C2*C3*C4*C5^2*S1*S4*nx^2 + 2*C1*C2*C3*C4*C5^2*S1*S4*ny^2 + 4*C1^2*C2*C3*C4*C5^2*S4*nx*ny + 2*C1*C4*C5^2*S1*S2*S3*S4*nx^2 - 2*C1*C4*C5^2*S1*S2*S3*S4*ny^2 - 8*C2^2*C3^2*C4*C5*S1*S5*nx*ny - 4*C1^2*C4*C5^2*S2*S3*S4*nx*ny - 4*C1^2*C2*C3^2*C4*C5*S2*S5*nx^2 - 4*C1^2*C2^2*C3*C4*C5*S3*S5*nx^2 + 4*C1^2*C2*C3^2*C4*C5*S2*S5*ny^2 + 4*C1^2*C2^2*C3*C4*C5*S3*S5*ny^2 - 2*C1*C2*C3*C5*S4*S5*nx*ny + 4*C1*C2*C4*C5*S1*S2*S5*nx*ny + 4*C1*C3*C4*C5*S1*S3*S5*nx*ny + 8*C1*C2*C3*C4*C5*S2*S3*S5*nx^2 + 4*C1*C2*C3*C4^2*S1*S2*S3*nx*ny - 4*C1*C2*C3*C5^2*S1*S2*S3*nx*ny - 4*C1*C2*C3*C4^2*C5^2*S1*S2*S3*nx*ny + 8*C2*C3*C4*C5*S1*S2*S3*S5*nx*ny - 8*C1*C2*C3^2*C4*C5*S1*S2*S5*nx*ny - 8*C1*C2^2*C3*C4*C5*S1*S3*S5*nx*ny - 1)^(1/2)) - atan2(C1*C5*S4*ny - C2*C3*S5*nx - C5*S1*S4*nx + S2*S3*S5*nx - C2*S1*S3*S5*ny - C3*S1*S2*S5*ny - C2*C4*C5*S3*nx - C3*C4*C5*S2*nx - C1*C2*S3*S5*nx - C1*C3*S2*S5*nx + C1*C2*C3*C4*C5*nx + C2*C3*C4*C5*S1*ny - C1*C4*C5*S2*S3*nx - C4*C5*S1*S2*S3*ny, C1*C4*ny - C4*S1*nx + C2*S3*S4*nx + C3*S2*S4*nx + C1*S2*S3*S4*nx + S1*S2*S3*S4*ny - C1*C2*C3*S4*nx - C2*C3*S1*S4*ny);
%Solution 2
else
    q1 = atan2(ty, tx);
    S1 = sin(q1);
    C1 = cos(q1);
 
    q2 = atan2(- C1^2*tx^2 + C1^2*ty^2 - 2*S1*C1*tx*ty - L1^2 + L2^2 - ty^2 - tz^2, (4*C1^6*L1^2*ty^2 - 12*C1^4*L1^2*ty^2 - C1^4*tx^4 + 6*C1^4*tx^2*ty^2 - C1^4*ty^4 + 8*C1^3*L1^2*S1*tx*ty - 4*C1^3*S1*tx^3*ty + 4*C1^3*S1*tx*ty^3 + 2*C1^2*L1^2*tx^2 + 10*C1^2*L1^2*ty^2 + 2*C1^2*L2^2*tx^2 - 2*C1^2*L2^2*ty^2 - 6*C1^2*tx^2*ty^2 - 2*C1^2*tx^2*tz^2 + 2*C1^2*ty^4 + 2*C1^2*ty^2*tz^2 + 8*C1*L1^2*S1^3*tx*ty - 4*C1*L1^2*S1*tx*ty + 4*C1*L2^2*S1*tx*ty - 4*C1*S1*tx*ty^3 - 4*C1*S1*tx*ty*tz^2 - L1^4 + 2*L1^2*L2^2 + 4*L1^2*S1^6*ty^2 - 2*L1^2*ty^2 + 2*L1^2*tz^2 - L2^4 + 2*L2^2*ty^2 + 2*L2^2*tz^2 - ty^4 - 2*ty^2*tz^2 - tz^4)^(1/2)) - angle(-L1*(ty*C1^2*S1*i + tx*C1*i + ty*S1^3*i - tz));
    S2 = sin(q2);
    C2 = cos(q2);
 
    q3 = atan2(C1*C2*tx - S2*tz - L1 + C2*S1*ty, C2*tz + C1*S2*tx + S1*S2*ty);
    S3 = sin(q3);
    C3 = cos(q3);
 
    q4 = atan2(C1*ay - S1*ax, C1*C2*C3*ax - C3*S2*az - C2*S3*az + C2*C3*S1*ay - C1*S2*S3*ax - S1*S2*S3*ay);
    S4 = sin(q4);
    C4 = cos(q4);
 
    q5 = atan2(C2*C3*nx - S2*S3*nx + C1*C2*S3*nx + C1*C3*S2*nx + C2*S1*S3*ny + C3*S1*S2*ny, S1*S4*nx - C1*S4*ny + C2*C4*S3*nx + C3*C4*S2*nx + C4*S1*S2*S3*ny - C1*C2*C3*C4*nx - C2*C3*C4*S1*ny + C1*C4*S2*S3*nx);
    S5 = sin(q5);
    C5 = cos(q5);
 
    q6 = atan2(1, (nx^2 + ny^2 + C1^2*nx^2 + C4^2*nx^2 - C1^2*ny^2 - C4^2*ny^2 - 2*C1^2*C4^2*nx^2 - C1^2*C5^2*nx^2 - C2^2*C4^2*nx^2 + C2^2*C5^2*nx^2 - C3^2*C4^2*nx^2 + C3^2*C5^2*nx^2 - C4^2*C5^2*nx^2 + 2*C1^2*C4^2*ny^2 + C1^2*C5^2*ny^2 + C2^2*C4^2*ny^2 - C2^2*C5^2*ny^2 + C3^2*C4^2*ny^2 - C3^2*C5^2*ny^2 + C4^2*C5^2*ny^2 + C1^2*C2^2*C4^2*nx^2 - C1^2*C2^2*C5^2*nx^2 + C1^2*C3^2*C4^2*nx^2 - C1^2*C3^2*C5^2*nx^2 + 2*C2^2*C3^2*C4^2*nx^2 + 2*C1^2*C4^2*C5^2*nx^2 - 2*C2^2*C3^2*C5^2*nx^2 + C2^2*C4^2*C5^2*nx^2 + C3^2*C4^2*C5^2*nx^2 - C1^2*C2^2*C4^2*ny^2 + C1^2*C2^2*C5^2*ny^2 - C1^2*C3^2*C4^2*ny^2 + C1^2*C3^2*C5^2*ny^2 - 2*C2^2*C3^2*C4^2*ny^2 - 2*C1^2*C4^2*C5^2*ny^2 + 2*C2^2*C3^2*C5^2*ny^2 - C2^2*C4^2*C5^2*ny^2 - C3^2*C4^2*C5^2*ny^2 + 2*C1*S1*nx*ny - 2*C1*C4*C5*S5*nx^2 - 4*C1*C4^2*S1*nx*ny - 2*C1*C5^2*S1*nx*ny - 2*C1*C2*C4^2*S2*nx^2 + 2*C1*C2*C5^2*S2*nx^2 - 2*C1*C3*C4^2*S3*nx^2 + 2*C1*C3*C5^2*S3*nx^2 - 2*C1^2*C2^2*C3^2*C4^2*nx^2 + 2*C1^2*C2^2*C3^2*C5^2*nx^2 - C1^2*C2^2*C4^2*C5^2*nx^2 - C1^2*C3^2*C4^2*C5^2*nx^2 - 2*C2^2*C3^2*C4^2*C5^2*nx^2 + 2*C1^2*C2^2*C3^2*C4^2*ny^2 - 2*C1^2*C2^2*C3^2*C5^2*ny^2 + C1^2*C2^2*C4^2*C5^2*ny^2 + C1^2*C3^2*C4^2*C5^2*ny^2 + 2*C2^2*C3^2*C4^2*C5^2*ny^2 + 2*C1*C2^2*C4^2*S1*nx*ny - 2*C1*C2^2*C5^2*S1*nx*ny + 2*C1*C3^2*C4^2*S1*nx*ny - 2*C1*C3^2*C5^2*S1*nx*ny + 4*C1*C4^2*C5^2*S1*nx*ny + 4*C1*C2*C3^2*C4^2*S2*nx^2 - 4*C1*C2*C3^2*C5^2*S2*nx^2 + 4*C1*C2^2*C3*C4^2*S3*nx^2 + 2*C1*C2*C4^2*C5^2*S2*nx^2 - 4*C1*C2^2*C3*C5^2*S3*nx^2 + 2*C1*C3*C4^2*C5^2*S3*nx^2 + 2*C2*C3*C4*S4*nx*ny + 2*C1^2*C2^2*C3^2*C4^2*C5^2*nx^2 - 2*C1^2*C2^2*C3^2*C4^2*C5^2*ny^2 - 2*C4*C5*S1*S5*nx*ny - 2*C4*S2*S3*S4*nx*ny - 2*C2*C4*C5*S2*S5*nx^2 - 2*C3*C4*C5*S3*S5*nx^2 + 2*C2*C4*C5*S2*S5*ny^2 + 2*C3*C4*C5*S3*S5*ny^2 - 2*C2*C4*S1*S3*S4*nx^2 - 2*C3*C4*S1*S2*S4*nx^2 - 2*C2*C4^2*S1*S2*nx*ny + 2*C2*C5^2*S1*S2*nx*ny - 2*C3*C4^2*S1*S3*nx*ny + 2*C3*C5^2*S1*S3*nx*ny + 4*C1*C2^2*C4*C5*S5*nx^2 + 4*C1*C3^2*C4*C5*S5*nx^2 - 2*C2*C3*C4^2*S2*S3*nx^2 + 2*C2*C3*C5^2*S2*S3*nx^2 + 2*C2*C3*C4^2*S2*S3*ny^2 - 2*C2*C3*C5^2*S2*S3*ny^2 + 4*C2^2*C4*C5*S1*S5*nx*ny + 4*C3^2*C4*C5*S1*S5*nx*ny - 2*C5*S1*S2*S3*S4*S5*nx^2 + 4*C1^2*C4*S2*S3*S4*nx*ny + 2*C4*C5^2*S2*S3*S4*nx*ny + 2*C1^2*C2*C4*C5*S2*S5*nx^2 + 4*C2*C3^2*C4*C5*S2*S5*nx^2 + 2*C1^2*C3*C4*C5*S3*S5*nx^2 + 4*C2^2*C3*C4*C5*S3*S5*nx^2 - 2*C1^2*C2*C4*C5*S2*S5*ny^2 - 4*C2*C3^2*C4*C5*S2*S5*ny^2 - 2*C1^2*C3*C4*C5*S3*S5*ny^2 - 4*C2^2*C3*C4*C5*S3*S5*ny^2 + 2*C2*C4*C5^2*S1*S3*S4*nx^2 + 2*C3*C4*C5^2*S1*S2*S4*nx^2 + 4*C2*C3^2*C4^2*S1*S2*nx*ny - 4*C2*C3^2*C5^2*S1*S2*nx*ny + 4*C2^2*C3*C4^2*S1*S3*nx*ny + 2*C2*C4^2*C5^2*S1*S2*nx*ny - 4*C2^2*C3*C5^2*S1*S3*nx*ny + 2*C3*C4^2*C5^2*S1*S3*nx*ny - 8*C1*C2^2*C3^2*C4*C5*S5*nx^2 + 2*C1^2*C2*C3*C4^2*S2*S3*nx^2 - 2*C1^2*C2*C3*C5^2*S2*S3*nx^2 + 2*C2*C3*C4^2*C5^2*S2*S3*nx^2 - 2*C1^2*C2*C3*C4^2*S2*S3*ny^2 + 2*C1^2*C2*C3*C5^2*S2*S3*ny^2 - 2*C2*C3*C4^2*C5^2*S2*S3*ny^2 + 2*C1*C2*C4*S3*S4*nx*ny + 2*C1*C3*C4*S2*S4*nx*ny - 4*C1*C2^2*C3^2*C4^2*S1*nx*ny + 4*C1*C2^2*C3^2*C5^2*S1*nx*ny - 2*C1*C2^2*C4^2*C5^2*S1*nx*ny - 2*C1*C3^2*C4^2*C5^2*S1*nx*ny + 2*C2*C5*S3*S4*S5*nx*ny + 2*C3*C5*S2*S4*S5*nx*ny + 2*C1*C2*C3*C4*S1*S4*nx^2 - 2*C1*C2*C3*C4*S1*S4*ny^2 - 4*C1*C2*C3^2*C4^2*C5^2*S2*nx^2 - 4*C1*C2^2*C3*C4^2*C5^2*S3*nx^2 + 2*C2*C3*C5*S1*S4*S5*nx^2 - 4*C1^2*C2*C3*C4*S4*nx*ny - 2*C2*C3*C4*C5^2*S4*nx*ny - 2*C1*C4*S1*S2*S3*S4*nx^2 + 2*C1*C4*S1*S2*S3*S4*ny^2 - 4*C2*C3^2*C4^2*C5^2*S1*S2*nx*ny - 4*C2^2*C3*C4^2*C5^2*S1*S3*nx*ny + 2*C1*C5*S2*S3*S4*S5*nx*ny - 2*C1^2*C2*C3*C4^2*C5^2*S2*S3*nx^2 + 2*C1^2*C2*C3*C4^2*C5^2*S2*S3*ny^2 + 2*C1*C2*C5*S1*S3*S4*S5*nx^2 + 2*C1*C3*C5*S1*S2*S4*S5*nx^2 - 2*C1*C2*C5*S1*S3*S4*S5*ny^2 - 2*C1*C3*C5*S1*S2*S4*S5*ny^2 - 2*C1*C2*C4*C5^2*S3*S4*nx*ny - 2*C1*C3*C4*C5^2*S2*S4*nx*ny + 4*C1*C2^2*C3^2*C4^2*C5^2*S1*nx*ny - 4*C1^2*C2*C5*S3*S4*S5*nx*ny - 4*C1^2*C3*C5*S2*S4*S5*nx*ny - 2*C1*C2*C3*C4*C5^2*S1*S4*nx^2 + 2*C1*C2*C3*C4*C5^2*S1*S4*ny^2 + 4*C1^2*C2*C3*C4*C5^2*S4*nx*ny + 2*C1*C4*C5^2*S1*S2*S3*S4*nx^2 - 2*C1*C4*C5^2*S1*S2*S3*S4*ny^2 - 8*C2^2*C3^2*C4*C5*S1*S5*nx*ny - 4*C1^2*C4*C5^2*S2*S3*S4*nx*ny - 4*C1^2*C2*C3^2*C4*C5*S2*S5*nx^2 - 4*C1^2*C2^2*C3*C4*C5*S3*S5*nx^2 + 4*C1^2*C2*C3^2*C4*C5*S2*S5*ny^2 + 4*C1^2*C2^2*C3*C4*C5*S3*S5*ny^2 - 2*C1*C2*C3*C5*S4*S5*nx*ny + 4*C1*C2*C4*C5*S1*S2*S5*nx*ny + 4*C1*C3*C4*C5*S1*S3*S5*nx*ny + 8*C1*C2*C3*C4*C5*S2*S3*S5*nx^2 + 4*C1*C2*C3*C4^2*S1*S2*S3*nx*ny - 4*C1*C2*C3*C5^2*S1*S2*S3*nx*ny - 4*C1*C2*C3*C4^2*C5^2*S1*S2*S3*nx*ny + 8*C2*C3*C4*C5*S1*S2*S3*S5*nx*ny - 8*C1*C2*C3^2*C4*C5*S1*S2*S5*nx*ny - 8*C1*C2^2*C3*C4*C5*S1*S3*S5*nx*ny - 1)^(1/2)) - atan2(C1*C5*S4*ny - C2*C3*S5*nx - C5*S1*S4*nx + S2*S3*S5*nx - C2*S1*S3*S5*ny - C3*S1*S2*S5*ny - C2*C4*C5*S3*nx - C3*C4*C5*S2*nx - C1*C2*S3*S5*nx - C1*C3*S2*S5*nx + C1*C2*C3*C4*C5*nx + C2*C3*C4*C5*S1*ny - C1*C4*C5*S2*S3*nx - C4*C5*S1*S2*S3*ny, C1*C4*ny - C4*S1*nx + C2*S3*S4*nx + C3*S2*S4*nx + C1*S2*S3*S4*nx + S1*S2*S3*S4*ny - C1*C2*C3*S4*nx - C2*C3*S1*S4*ny);
    
end

q_vec = [q1,q2,q3,q4,q5,q6];

end

