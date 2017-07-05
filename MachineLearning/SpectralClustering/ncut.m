function [V1,D1]=ncut(W);
% [V1,D1]=ncut(W,nv);
%
% solve generalized eigenproblem Wy={\mu}My
%
% optional third argument specifies minimum allowed weight (e.g. 0.01)


m=sum(W,1);
N=length(m);
M=sparse(1:N,1:N,m);
B=inv(sqrt(M));

% solve generalized eigensystem
OPTIONS.tol=1e-4;
OPTIONS.maxit=20;
%OPTIONS.disp=0;
C=B*W*B;
format long
[V1,D1]=eig(C);

D1=1-D1;
V1=B*V1;

% sort the eigenvalues
[dvalues,dindex]=sort(diag(D1));
%keyboard;
V1=V1(:,dindex);
D1=diag(dvalues);
format
