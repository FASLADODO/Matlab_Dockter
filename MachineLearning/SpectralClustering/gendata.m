function data = gendata(mean,sigmax,sigmay,N);
data=repmat(mean,1,N)+twonorm(sigmax,sigmay,N);
