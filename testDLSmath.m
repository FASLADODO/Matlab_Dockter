n = 3;
m = 5;

x1 = randn(m,n);
x2 = randn(m,n);

y1 = randn(m,1);
y2 = randn(m,1);

b1 = randn(n,1);
b2 = randn(n,1);


2*x1'*y1
2*x1'*x1*b1

out1 = 2*x2'*x2*b1*(b2'*x2'*x2*b2)
out2 = 2*x2'*x2*b2*(b1'*x2'*x2*b1)


out3 = 2*x1'*y1*(2*b1'*x2'*y2)
out4 = 2*x2'*y2*(2*b1'*x1'*y1)