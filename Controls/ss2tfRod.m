%Tf to SS

%SS matrices

A = [0, 1, 0;
    0, 0, 1;
    -3, -4, -2];

B = [0;
    0;
    1];

C = [5, 1, 0];

D = 0;

syms s

TF = C*inv(s*eye(length(A)) - A)*B + D

