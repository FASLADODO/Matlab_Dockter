function data = twonorm(sigmax,sigmay,N)
X=normrnd(0,sigmax,1,N);
Y=normrnd(0,sigmay,1,N);
data=cat(1,X,Y);
