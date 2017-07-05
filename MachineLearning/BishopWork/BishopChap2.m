%%  Bishop chap 2

%Parametric distrubutions defined by small number of parameters

%binary variables

heads = 0.4;
tails = 0.6;
nn = 100;
rng('shuffle');  % for reproducibility
r = binornd(ones(nn,1),heads);

p1 = sum(r)/length(r)

bar([1,2],[p1,1-p1])

%% Binomial

N = 20;
mu = 0.25;

M = 0:6;

for m = M
   bin(m+1) = (factorial(N)/(factorial(N-m)*factorial(m)) )*(mu^m) *((1-mu)^(N-m));
end

plot(bin)

E = N*mu
varE = N*mu*(1-mu)
var(bin)

%% beta function

N = 20;
MU = 0:0.01:1;


a = 0.1;
b = 0.1;

for ii = 1:length(MU)
   T_ab = gamma(a+b)
   T_a = gamma(a );
   T_b = gamma(b );

   beta(ii) = (T_ab/(T_a+T_b))*(MU(ii)^(a-1)) *((1-MU(ii))^(b-1));
end

plot(MU,beta)

%% beta binomial prior

ph = 0.8;

N = 10;
MU = 0:0.01:1;
m = round ( ph*N );
l = N - m;

a = 2;
b = 2;

for ii = 1:length(MU)
   T_abml = gamma(a+b+l+m)
   T_ma = gamma(a +m);
   T_lb = gamma(b +l);

   beta(ii) = (T_abml/(T_ma+T_lb))*(MU(ii)^(m+a-1)) *((1-MU(ii))^(l+b-1));
end

p_x1 = (m+a)/(m+a+l+b)

plot(MU,beta)

%% Multinomial distributions 

% want to find u_k probabilities for each outcome

%x = [ 0,0,0,1,0,0] 1-of-k scheme
N=50;
k = 6;

x = oneOfK(N,k);

mk = sum(x,1) %number of observations

%maximum likelihood solutions
u_k = -mk/N

plot([1:k],u_k)


%% conditional gaussians (IT WORKS)
nn = 500;

mu1 = [2,3];
sigma1 = [2,0.5;0.5,2];
mu2 = [4,5];
sigma2 = [2,0.5;0.5,2];

%Regular distrubutions
x1 = mvnrnd(mu1,sigma1,nn);
x2 = mvnrnd(mu2,sigma2,nn);

mu1 = mean(x1);
sig1 = cov(x1)
mu2 = mean(x2);
sig2 = cov(x2)

p1 = gaussianProbMVarray(x1, sig1, mu1);
p2 = gaussianProbMVarray(x2, sig2, mu2);

%plot seperate
figure
scatter(x1(:,1),x1(:,2),'r.');
hold on
handle1 = Surface3D(x1(:,1),x1(:,2),p1);
hold on
scatter(x2(:,1),x2(:,2),'b.');
hold on
handle2 = Surface3D(x2(:,1),x2(:,2),p2);
hold off

[mua,mub,sigabbbinv,siga_b] = ConditionalGaussian(x1,x2);
[pa_b] = ConditionalGaussianDist(x2,mua,mub,sigabbbinv,siga_b);

% [cparamLG,csigL,csigin,css] = LinearGaussianTrain(x1,x2);
% [pa_b,cmeany] = LinearGaussianOnline(x1,x2,cparamLG,csigL,csigin,css);

%plot seperate
figure
scatter(x1(:,1),x1(:,2),'r.');
hold on
scatter(x2(:,1),x2(:,2),'b.');
hold on
handle2 = Surface3D(x2(:,1),x2(:,2),pa_b);
hold off

%% linear gaussian using conditionals Page 91 Bishop

noiz = 3;
params = [3;2.5];
range = [1,25];
nn = 500;

xo = linspace(range(1),range(2),nn)';
data = [ones(nn,1), xo];
y = data*params + randn(nn,1)*noiz;

figure
scatter(xo,y)

[mua,mub,sigabbbinv,siga_b] = ConditionalGaussian(data,y);
[PL] = ConditionalGaussianDist(data,mua,mub,sigabbbinv,siga_b);

%Train params
% [paramLG,sigL,sigin,ss] = LinearGaussianTrain(data,y);
% 
% [PL,meany] = LinearGaussianOnline(data,y,paramLG,sigL,sigin,ss);

figure
scatter(xo,y,'b.')
hold on
% plot(xo,meany,'r-')
% hold on
handle2 = Surface3D(xo,y,PL);
hold off


% %Now do it for grid
% limz = [min(xo),min(y); max(xo),max(y)];
% xlin = linspace(limz(1,1),limz(2,1),50)';
% ylin = linspace(limz(1,2),limz(2,2),50)';
% [XG,YG] = meshgrid(xlin,ylin);
% data2 = [ones(length(XG(:)),1), XG(:)];
% Y2 = YG(:);
% 
% [PLG,meanyg] = LinearGaussianOnline(data2,Y2,paramLG,sigL,sigin,ss);
% 
% figure
% scatter(xo,y,'b.')
% hold on
% plot(xo,meany,'r-')
% hold on
% handle3 = Surface3D(XG(:),Y2,PLG);
% hold off

%% linear gaussian for cube functions

noiz = 30;
params = [3;2.5;-0.5;-0.1];
range = [1,25];
nn = 500;

xo = linspace(range(1),range(2),nn)';
data = [ones(nn,1), xo, xo.^2, xo.^3];
y = data*params + randn(nn,1)*noiz;

figure
scatter(xo,y)


[PL,meany] = LinearGaussian(data,y,data,y);

figure
scatter(xo,y,'b.')
hold on
plot(xo,meany,'r-')
hold on
handle2 = Surface3D(xo,y,PL);
hold off

%Now do it for grid
limz = [min(xo),min(y); max(xo),max(y)];
xlin = linspace(limz(1,1),limz(2,1),50)';
ylin = linspace(limz(1,2),limz(2,2),50)';
[XG,YG] = meshgrid(xlin,ylin);
data2 = [ones(length(XG(:)),1), XG(:), XG(:).^2, XG(:).^3];
Y2 = YG(:);

[PLG,meanyg] = LinearGaussian(data,y,data2,Y2);

figure
scatter(xo,y,'b.')
hold on
plot(xo,meany,'r-')
hold on
handle3 = Surface3D(XG(:),Y2,PLG);
hold off


%% 3rd and 4th moments



nn = 100;

if(~exist('xtemp','var') )
    disp('new rand')
    xtemp = randn(nn,2);
end


mu = [4,4];
scale = [3,3];
x0 = [xtemp.*repmat(scale,nn,1) + repmat(mu,nn,1) ];

x1 = [x0(:,1), x0(:,1)*2.5 + randn(nn,1)*2]; %linear

%plot seperate
figure
scatter(x1(:,1),x1(:,2),'r.');
% hold on
% scatter(x2(:,1),x2(:,2),'b.');
hold off

%Mean
m1 = moment(x1,1)

%variance
m2 = moment(x1,2)

%Skewness
m3 = moment(x1,3)

%Kurtosis
m4 = moment(x1,4)

