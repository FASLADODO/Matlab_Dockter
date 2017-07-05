function [Abar,Bbar,Cbar] = KalmanCCF(A,B,C)
%%designed for A 4*4 and B 4*1
cont=B;
cont(1:4,2:2)=A*B;
cont(1:4,3:3)=A*A*B;
cont(1:4,4:4)=A*A*A*B;
[Q,S,V]=svd(cont);
Abar=inv(Q)*A*Q;
Bbar=inv(Q)*B;
Cbar=C*Q;