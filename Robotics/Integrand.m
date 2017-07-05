function K = Integrand()
%%Part C) program that creates exact K matrix for the plane element
syms x y
B=(1/24)*[-(4-y),0,4-y,0,4+y,0,-(4+y),0;0,-(1.5-x),0,-(1.5+x),0,(1.5+x),0,(1.5-x);-(1.5-x),-(4-y),-(1.5+x),4-y,1.5+x,4+y,1.5-x,-(4+y)];
D=(10^7/0.91)*[1,0.3,0;0.3,1,0;0,0,0.35];
Btran=transpose(B);
Z=Btran*D*B;
K1=int(Z,x);
K2=subs(K1,{'x'},{1.5})-subs(K1,{'x'},{-1.5});
K3=int(K2,y);
K=subs(K3,{'y'},{4})-subs(K3,{'y'},{-4});
