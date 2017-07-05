%Hoeffding Inequality

nn = 1000;

rng('shuffle')
huhin = randn(nn,1) + 0.8;

Zi = huhin > 0;

z0 = Zi == 0;
P0 = sum(z0)/length(z0);

z1 = Zi == 1;
P1 = sum(z1)/length(z1);


phi = P1



rng('shuffle')
huhout = randn(nn,1) + 0.8;

scatter(1:length(huhout),huhout,'b.')

acc = huhout > 0;

phihat = mean(acc)


gamma = 0.01
ineq1 = abs(phi - phihat)
ineq2 = 2*exp(-2*(gamma^2)*nn)
