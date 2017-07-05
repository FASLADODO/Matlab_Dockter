Data = randn(100,10);
X_World = randn(100,1);
Y_World = randn(100,1);

[BIC_Results] = BICRMSE(Data,X_World,Y_World);

figure
plot(BIC_Results.BestRMSE)