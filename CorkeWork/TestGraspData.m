%loadings
load test6.mat

%params
Phi1 = [-8.02086003501975e-08;
1.55286905458007e-05;
0.00795344249784777;
8.23733045077119;
-0.00299236282304311];

Phi2 = [-1.24436831310503e-08;
1.23673348605010e-05;
0.00652545188345528;
6.75893262890734;
-0.00228098997419065];

classify = [];
est_class = 0;
cntr = 0;
errors{1} = [];

Dx = [science(:,4),science(:,3),science(:,2),science(:,2).^2,science(:,2).^3];

e1 = abs(science(:,5) - Dx*Phi1);
e2 = abs(science(:,5) - Dx*Phi2);

if(sum(e1) < sum(e2))
    est_class = 1;
else
    est_class = 2;
end

cummy = [cumsum(e1),cumsum(e2)];

figure(2)
plot(e1,'rx')
hold on
plot(science(:,7),'bo')
hold off

title('science 1')

figure(4)
plot(e2,'rx')
hold on
plot(science(:,8),'bo')
hold off

title('science 2')

figure(5)
plot(science(:,5),'rx')


title('Forces')