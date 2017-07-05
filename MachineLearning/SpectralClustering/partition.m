W=[ 0  20 1 1 ; 
    20 0 1  1; 
    1 1 0 20 ;
    1 1 20 0;];

% lets create a 50x50 graph
% first 25 vertices are on one side
% the other 25 are on the other side.

A= rand(25)/4.0+0.5;
B=  randn(25)/4.0+ 0.5;
C=0.5-randn(25)/4.0;

A = ones(25);
B= ones(25);
C=zeros(25);
W=[ A C;C' B]'

L=W-diag(sum(W));
[V,D]=eig(L);
