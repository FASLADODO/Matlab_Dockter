
xt= -10:0.1:10;

figure 
scatter(xt,erf(xt),'r.')
title('this is how erf works')

%So erf is always 1 for x > 1


%% 1D integral

% make some data
nn = 100;
a = randn(nn,1);

%make the function
func = @(x) sum(exp(- norm(x-a) ) );

%some bounds
xmin = -5;
xmax = 5;
xtemp = -5:0.1:5;

%compute intergal
intval = integral(func,xmin,xmax,'ArrayValued',true);

%get rbf at each point
for ii = 1:length(xtemp);
    prob(ii) = func(xtemp(ii));
end

%make plot
figure
scatter(a,zeros(nn,1),'ro')
hold on
plot(xtemp,prob,'b')
hold off
str = sprintf('1D data with integral = %d', intval);
title(str)

%% 2D integral

a = randn(10,2);

func = @(x) exp(- sum(pdist2(x,a),2) ) ;

func(randn(10,2))
%% try the func integral idea

nn = 100;

D1 = randn(nn,2);
D2 = randn(nn,2) + repmat([3,3],nn,1);

figure
scatter(D1(:,1),D1(:,2),'r.')
hold on
scatter(D2(:,1),D2(:,2),'b.')
hold off

Dwithin = D1;
Dbetween = D2;

DAll = [Dwithin;Dbetween];

Bounds = DataBounds(DAll)*3;

%scale by lengths
nw = size(Dwithin,1);
nb = size(Dbetween,1);

xmin = Bounds(1,:);
xmax = Bounds(2,:);


normmin = NormRowWise(Dwithin - repmat(xmin,nw,1));
normmax = NormRowWise(Dwithin - repmat(xmax,nw,1));

%Wolfram alpha: int(exp(-norm(x)))
radbasint = @(x) 0.5.*exp(-x).*( -(exp(x) - 1).^2 + 2.*exp(x) + exp(2*x) - 1);
radbasint(5) - radbasint(2)

intgrl = sum(radbasint(normmax) - radbasint(normmin))


%% Try the average of dx*dy space THIS IS THE GOOD ONE!!!!!!!!!!!!!!!!
% 1D

nn = 2000;
S = 1;
mu = [4];

Dwithin = randn(nn,S)*2;


%scale by lengths
[non,son] = size(Dwithin);


nw = non;
sigw = norm(std(Dwithin));
bandwidthw = 1.06*sigw*(nw^(-1/5)) ;

%get pdists
distwithin = pdist2(Dwithin,Dwithin);

%get rbfs
rw = exp( -(distwithin.^2)./(2*bandwidthw^2) );

%get within and between class unsca;ed
Pwithin_us = sum(rw,2); %/(noff^2);%

%integral estimate
scalew = ScaleRBF(Dwithin,Pwithin_us); %prod(stepsize_within)*sum(Pwithin_us)

%scale the probs
Pwithin = Pwithin_us/scalew;


%compare with mvnpdf
MU = mean(Dwithin);
SIGMA = cov(Dwithin);
PM = mvnpdf(Dwithin,MU,SIGMA);


figure
scatter(Dwithin(:,1),zeros(nn,1),'ro')
hold on
scatter(Dwithin(:,1),PM,'g.')
hold off
title('Matlabs mvnpdf')


figure
scatter(Dwithin(:,1),zeros(nn,1),'ro')
hold on
scatter(Dwithin(:,1),Pwithin,'g.')
hold off
trapw = trapz(Dwithin,Pwithin_us)
title('Rods PDF')


ScaleRBF(Dwithin,Pwithin)


%% 2D version

nn = 2000;
S = 2;
mu = [3,3];

Dwithin = randn(nn,S)*2;

%scale by lengths
[non,son] = size(Dwithin);

nw = non;
sigw = norm(std(Dwithin));
bandwidthw = 1.06*sigw*(nw^(-1/5)) ;

%get pdists
distwithin = pdist2(Dwithin,Dwithin);

%get rbfs
rw = exp( -(distwithin.^2)./(2*bandwidthw^2) );

%get within and between class unsca;ed
Pwithin_us = sum(rw,2); %/(noff^2);%

[scale2,Grid,p] = integralRBF(Dwithin,bw);
scalew = ScaleRBF(Dwithin,Pwithin_us); %prod(stepsize_within)*sum(Pwithin_us)

%compare scales
scalew
scale2

%scale the probs
Pwithin = Pwithin_us/scalew;

%compare with mvnpdf
MU = mean(Dwithin);
SIGMA = cov(Dwithin);
PM = mvnpdf(Dwithin,MU,SIGMA);


figure
scatter(Dwithin(:,1),Dwithin(:,2),'ro')
hold on
Surface3D(Dwithin(:,1),Dwithin(:,2),PM);
hold off
title('Matlabs mvnpdf')


figure
scatter(Dwithin(:,1),Dwithin(:,2),'ro')
hold on
Surface3D(Dwithin(:,1),Dwithin(:,2),Pwithin);
hold off
title('Rods PDF')

figure
scatter(Grid(:,1),Grid(:,2),'ro')
hold on
Surface3D(Grid(:,1),Grid(:,2),p);
hold off
title('grid PDF us')

figure
scatter(Grid(:,1),Grid(:,2),'ro')
hold on
Surface3D(Grid(:,1),Grid(:,2),p/scale2);
hold off
title('grid PDF scaled')

ScaleRBF(Dwithin,Pwithin)

%% 3D version

nn = 2000;
S = 3;
mu = [3,3,3];

Dwithin = randn(nn,S)*2;

%scale by lengths
[non,son] = size(Dwithin);

nw = non;
sigw = norm(std(Dwithin));
bandwidthw = 1.06*sigw*(nw^(-1/5)) ;

%get pdists
distwithin = pdist2(Dwithin,Dwithin);

%get rbfs
rw = exp( -(distwithin.^2)./(2*bandwidthw^2) );

%get within and between class unsca;ed
Pwithin_us = sum(rw,2); %/(noff^2);%

scalew = ScaleRBF(Dwithin,Pwithin_us); %prod(stepsize_within)*sum(Pwithin_us)

%scale the probs
Pwithin = Pwithin_us/scalew;

%compare with mvnpdf
MU = mean(Dwithin);
SIGMA = cov(Dwithin);
PM = mvnpdf(Dwithin,MU,SIGMA);


figure
scatter3(Dwithin(:,1),Dwithin(:,2),Dwithin(:,3),10,PM)
colormap cool
colorbar
title('Matlabs mvnpdf')


figure
scatter3(Dwithin(:,1),Dwithin(:,2),Dwithin(:,3),10,Pwithin)
colormap cool
colorbar
title('Rods pdf')

ScaleRBF(Dwithin,Pwithin)


%MAX PERCENTAGE SCALES EXACTLY AS 10^(-(ss+1))

%% Try figuring out an ideal rule:

nn = 2000;

storecoeffs = [];
svals = 1:9

%run through number of dimensions
for S = svals

    storallthescales = [];
    storeit = [];
    storesig = [];
    scalez = [1,1.1,1.5,1.8,2,2.5,3,4,5,6];

    %run through number of scales
    for ii = 1:length(scalez)

        Dwithin = randn(nn,S)*scalez(ii);
        Dbetween = randn(nn,S);

        %scale by lengths
        [non,son] = size(Dwithin);
        [noff,soff] = size(Dbetween);

        nw = non;
        sigw = norm(std(Dwithin));
        bandwidthw = 1.06*sigw*(nw^(-1/5)) ;
        storesig = [storesig; sigw/sqrt(S)];
        
        sigmaactual = sigw/sqrt(S);

        nb = noff;
        sigb = norm(std(Dbetween));
        bandwidthb = 1.06*sigb*(nb^(-1/5)) ;

        invs = (nb^(5))/(1.06*sigw);

        %get pdists
        distwithin = pdist2(Dwithin,Dwithin);
        distbetween = pdist2(Dbetween,Dbetween);

        %get rbfs
        rw = exp( -(distwithin.^2)./(2*bandwidthw^2) );
        rb = exp( -(distbetween.^2)./(2*bandwidthb^2) );
        
%         scalew = 1*(non^2)*exp(- 0.8*son);
%         scaleb = 1*(noff^2)*exp(- 0.8*soff);

        scalew = 0.5*(non^2)*son^(-2.5);
        scaleb = 0.5*(noff^2)*soff^(-2.5);
        %get within and between class
        Pwithin = sum(rw,2)/scalew; %/(noff^2);%
        Pbetween = sum(rb,2)/scaleb;

        sw = sum(Pwithin);
        sb = sum(Pbetween);
        storeit = [storeit;  scalez(ii), sw];

        
        storallthescales = [storallthescales; S, sw];
        
    end

%     figure
%     plot(storeit(:,1),storeit(:,2))
%     str = sprintf('sigma vs sum dim, = %d',S);
%     title(str);
%     xlabel('simga')
%     ylabel('sum')
    
%     figure
%     plot(scalez,storesig)
%     str = sprintf('sigma vs sigma est, dim = %d',S);
%     title(str);
%     xlabel('simga')
%     ylabel('sigmaest')

%     f = fit(storallthescales(:,1),storallthescales(:,2),'exp1');

    storecoeffs = [ storecoeffs; storallthescales];
end

% figure
% plot(svals,storecoeffs(:,1))
% title('dim vs coeff1');
% xlabel('dim')
% ylabel('scale')


% f = fit(svals',storecoeffs(:,1),'exp1')


figure
plot(storecoeffs(:,1),storecoeffs(:,2))
xlabel('dim')
ylabel('sum(p)')

nn^2/2

fiteq = 'a*(x^b)';
f = fit(storecoeffs(:,1),storecoeffs(:,2),fiteq)

%% try summing the overall integral for a grid

Data = [Dwithin;Dbetween];

%scale by lengths
nw = size(Dwithin,1);
nb = size(Dbetween,1);

bounds = DataBounds(Data)*2
Grid = ndimgrid(bounds,1000);

%compute optimal bandwidth if not supplied

nall = sum([nw,nb]);
sig = norm(std(Data));
bandwidth = 1.06*sig*(nall^(-1/5)) ;

gw = pdist2(Grid,Dwithin);
gb = pdist2(Grid,Dbetween);

grw = exp( -(gw.^2)./(2*bandwidth^2) );
grb = exp( -(gb.^2)./(2*bandwidth^2) );

scalew = sum(sum(grw,2))
scaleb = sum(sum(grb,2))

probw = sum(grw,2)/scalew;
probb = sum(grb,2)/scaleb;

figure
scatter(D1(:,1),D1(:,2),'r.')
hold on
scatter(D2(:,1),D2(:,2),'b.')
hold on
Surface3D(Grid(:,1),Grid(:,2),probw);
hold off


%% try summing the overall integral for just the known data set

Data = [Dwithin;Dbetween];

%scale by lengths
nw = size(Dwithin,1);
nb = size(Dbetween,1);

%compute optimal bandwidth if not supplied
nall = sum([nw,nb]);
sig = norm(std(Data));
bandwidth = 1.06*sig*(nall^(-1/5)) ;


%get pdists
distwithin = pdist2(Dwithin,Dwithin);
distbetween = pdist2(Dbetween,Dbetween);
%get rbfs
rw = exp( -(distwithin.^2)./(2*bandwidth^2) );
rb = exp( -(distbetween.^2)./(2*bandwidth^2) );
%get within and between class
scalew = sum(sum(rw,2));
scaleb = sum(sum(rb,2));

probw = sum(rw,2)/scalew;
probb = sum(rb,2)/scaleb;

figure
scatter(D1(:,1),D1(:,2),'r.')
hold on
scatter(D2(:,1),D2(:,2),'b.')
hold on
Surface3D(D1(:,1),D1(:,2),probw);
hold on
Surface3D(D2(:,1),D2(:,2),probb);
hold off


