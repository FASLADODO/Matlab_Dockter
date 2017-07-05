%Weighted least squares


%%
nn = 1000;
noiz = 0.6;

rangen = [2 40];

params1 = [2 0.7 3]; %m, b


xm = linspace(rangen(1),rangen(2),nn)';

yx1 = xm.*params1(1) + rand(nn,1)*noiz.*NormRowWise(xm);

X1 = [xm, yx1];
Y1 = xm*params1(2) + yx1*params1(3) + rand(nn,1)*noiz.*NormRowWise(xm);


figure
scatter3(X1(:,1),X1(:,2),Y1,'r.')

title('Sample Data','FontSize', 12)
xlabel('state1','FontSize', 12)
ylabel('state2','FontSize', 12)


%% Get info on the datars

mu1 = mean(X1)

sig1 = cov(X1)

LS_params = pinv(X1)*Y1

Resid = abs(Y1 - X1*LS_params);

figure
plot(Resid)

weights = 1./Resid;

WB = diag(weights);

[NN,SS] = size(X1);

tic
WLS_param = inv(X1'*inv(WB)*X1)*X1'*inv(WB)*Y1;
toc

params1([2,3])
LS_params
WLS_param

[bw,sew_b,msew] = lscov(X1,Y1,WB)
