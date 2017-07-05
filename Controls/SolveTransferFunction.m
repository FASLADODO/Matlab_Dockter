%Solve A*X = B

syms s F

%X = [X1;X2;X3];

A = [2,-2,0;
    -2,(5*s+2),-5*s;
    0,-5*s,(10*s^2+7*s)];

B = [F;0;0];

X = A\B