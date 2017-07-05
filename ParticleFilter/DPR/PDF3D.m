%test continuous distributions

%http://www.mathworks.com/help/stats/fitgmdist.html

%% 2D

nn = 1000;
az = 20;
el = 30;


MU1 = [2 2];
SIGMA1 = [2 0; 0 .5];

MU2 = [-3 -3];
SIGMA2 = [1 0; 0 3;];

MU3 = [5 0];
SIGMA3 = [2 0; 0 1;];

X1 = mvnrnd(MU1,SIGMA1,nn);
X2 = mvnrnd(MU2,SIGMA2,nn);
X3 = mvnrnd(MU3,SIGMA3,nn);

X_All = [X1;X2;X3];
X12 = [X1;X2];
X13 = [X1;X3];
X23 = [X2;X3];

%get gaussian mixture model
GMModel1 = fitgmdist(X1,1);
GMModel2 = fitgmdist(X2,1);
GMModel3 = fitgmdist(X3,1);
GMModel12 = fitgmdist(X12,1);
GMModel13 = fitgmdist(X13,1);
GMModel23 = fitgmdist(X23,1);

limz = [min(X_All(:,1)) max(X_All(:,1)) min(X_All(:,2)) max(X_All(:,2))];
NP = 50;
xo1 = linspace(limz(1),limz(2),NP);
xo2 = linspace(limz(3),limz(4),NP);
[XGt1,XGt2] = meshgrid(xo1,xo2);
XG1 = reshape(XGt1,[],1);
XG2 = reshape(XGt2,[],1);


%In column form
Xon = [XG1,XG2];

F1 = pdf(GMModel1,Xon);
F2 = pdf(GMModel2,Xon);
F3 = pdf(GMModel3,Xon);
F12 = pdf(GMModel12,Xon);
F13 = pdf(GMModel13,Xon);
F23 = pdf(GMModel23,Xon);

Scale1 = 1/max(F1);
Scale2 = 1/max(F2);
Scale3 = 1/max(F3);
Scale12 = 1/max(F12);
Scale13 = 1/max(F13);
Scale23 = 1/max(F23);

F1 = F1 .*Scale1;
F2 = F2 .*Scale2;
F3 = F3 .*Scale3;
F12 = F12 .*Scale12;
F13 = F13 .*Scale13;
F23 = F23 .*Scale23;

figure
p1 = scatter(X1(:,1),X1(:,2),2,'r.');
hold on
m1 = Surface3Dalt(Xon(:,1),Xon(:,2),F1,limz,NP);
set(m1,'facecolor','none')
hold on

p2 = scatter(X2(:,1),X2(:,2),2,'g.');
hold on
m2 = Surface3Dalt(Xon(:,1),Xon(:,2),F2,limz,NP);
set(m2,'facecolor','none')

p3 = scatter(X3(:,1),X3(:,2),2,'b.');
hold on
m3 = Surface3Dalt(Xon(:,1),Xon(:,2),F3,limz,NP);
set(m3,'facecolor','none')
hold off

colormap(cool);
colorbar;
title('Scatter Plot and Fitted PDF ')
xlabel('x1')
ylabel('x2')
zlabel('probs')
legend([p1(1), p2(1), p3(1)],'class1','class2','class3')
axis([limz 0 1 ])
view(az, el);


figure
p1 = scatter(X1(:,1),X1(:,2),2,'r.');
hold on
m1 = Surface3Dalt(Xon(:,1),Xon(:,2),F13,limz,NP);
set(m1,'facecolor','none')
hold on

p3 = scatter(X3(:,1),X3(:,2),2,'b.');

hold off

colormap(cool);
colorbar;
title('Combined PDF ')
xlabel('x1')
ylabel('x3')
zlabel('probs')
legend([p1(1), p3(1)],'class1','class2')
axis([limz 0 1 ])
view(az, el);

%%

%Try spatial kullback leibler divergence  (ie compute KL in a window
%surrounding each spot

%"In other words, it is the amount of information lost when Q is used to
%approximate P"
%KL divergence is how good of an estimate P is of Q, in other words how
%similar are the data sets
%A low KL divergence means two data sets are very similar, high KL means
%high seperability

%over KL
D_KL_12 = sum(F1 .* log(F1./F2) )
D_KL_13 = sum(F1 .* log(F1./F3) )
D_KL_23 = sum(F2 .* log(F2./F3) )


ws = 3;
pp = 10;

%%ratios of entropy
KL_Plot_12 = [];
KL_Plot_13 = [];
KL_Plot_23 = [];
KL_Plot_21 = [];
KL_Plot_31 = [];
KL_Plot_32 = [];

%Loop through limits of data and compute probabilities in each windows.
%Then compute KL in each
limit = limz;
for ii = limit(1):limit(2) %X axis
    for jj = limit(3):limit(4) %Y axis
        
        xl1 = linspace(ii-ws/2, ii+ws/2, pp);
        xl2 = linspace(jj-ws/2, jj+ws/2, pp);
        [xg1,xg2] = meshgrid(xl1,xl2);
        XGC1 = reshape(xg1,[],1);
        XGC2 = reshape(xg2,[],1);
        
        XWIN = [XGC1,XGC2];
        
        F1_temp =  pdf(GMModel1,XWIN);
        F2_temp =  pdf(GMModel2,XWIN);
        F3_temp =  pdf(GMModel3,XWIN);
        F12_temp = pdf(GMModel12,XWIN);
        F13_temp = pdf(GMModel13,XWIN);
        F23_temp = pdf(GMModel23,XWIN);

        F1_temp = F1_temp .*Scale1;
        F2_temp = F2_temp .*Scale2;
        F3_temp = F3_temp .*Scale3;
        F12_temp = F12_temp .*Scale12;
        F13_temp = F13_temp .*Scale13;
        F23_temp = F23_temp .*Scale23;
        
        
%         KL_12_t = sum(F1_temp .* log(F1_temp./F2_temp) )
%         KL_21_t = sum(F2_temp .* log(F2_temp./F1_temp) )
%         KL_13_t = sum(F1_temp .* log(F1_temp./F3_temp) )
%         KL_31_t = sum(F3_temp .* log(F3_temp./F1_temp) )
%         KL_23_t = sum(F2_temp .* log(F2_temp./F3_temp) )
%         KL_32_t = sum(F3_temp .* log(F3_temp./F2_temp) )
        KL_12_t = -sum(F1_temp .* log2(F1_temp) );
        KL_21_t = -sum(F2_temp .* log2(F2_temp) );
        KL_13_t = -sum(F1_temp .* log2(F1_temp) );
        KL_31_t = -sum(F3_temp .* log2(F3_temp) );
        KL_23_t = -sum(F2_temp .* log2(F2_temp) );
        KL_32_t = -sum(F3_temp .* log2(F3_temp) );
        
        KL_12_c = -sum(F12_temp .* log2(F12_temp) );
        KL_13_c = -sum(F13_temp .* log2(F13_temp) );
        KL_23_c = -sum(F23_temp .* log2(F23_temp) );
        
        ratio_12 = (KL_12_t+KL_21_t)/KL_12_c;
        ratio_13 = (KL_13_t+KL_31_t)/KL_13_c;
        ratio_23 = (KL_23_t+KL_32_t)/KL_23_c;

        KL_Plot_12 = [KL_Plot_12; ii, jj, ratio_12];
        KL_Plot_13 = [KL_Plot_13; ii, jj, ratio_13];
        KL_Plot_23 = [KL_Plot_23; ii, jj, ratio_23];
    end
end

figure
s1 = scatter(X1(:,1),X1(:,2),2,'b.');
hold on
s3 = scatter(X3(:,1),X3(:,2),2,'r.');
hold on
k13 = Surface3Dalt(KL_Plot_13(:,1),KL_Plot_13(:,2),KL_Plot_13(:,3),limit,NP);
set(k13,'facecolor','none')
hold off

colormap(cool);
colorbar;
title('Scatter Plot and KL 1-3 ')
xlabel('x1')
ylabel('x3')
zlabel('D_KL')
legend([s1(1), s3(1)],'class1','class2')
axis([limz 0 10 ])
view(az, el);

%% 3D

nn = 1000;

MU1 = [2 2 2];
SIGMA1 = [2 0 0; 0 .5 0; 0 0 1];

MU2 = [-3 -3 -3];
SIGMA2 = [1 0 0; 0 3 0; 0 0 2];

MU3 = [0 4 10];
SIGMA3 = [2 0 0; 0 1 0; 0 0 .5];

X1 = mvnrnd(MU1,SIGMA1,nn);
X2 = mvnrnd(MU2,SIGMA2,nn);
X3 = mvnrnd(MU3,SIGMA3,nn);

figure
scatter3(X1(:,1),X1(:,2),X1(:,3),10,'r.')
hold on
scatter3(X2(:,1),X2(:,2),X2(:,3),10,'b.')
hold on
scatter3(X3(:,1),X3(:,2),X3(:,3),10,'g.')
hold off

title('Training Distributions')
xlabel('x')
ylabel('y')
zlabel('z')
legend('Group1','Group2','Group3')

X_All = [X1;X2;X3];

%get gaussian mixture model
GMModel = fitgmdist(X_All,3);

F1 = pdf(GMModel,X_All);

figure
scatter3(X_All(:,1),X_All(:,2),X_All(:,3),8,F1(:));

colormap(cool);
colorbar;
title('Density Estimates')
xlabel('x')
ylabel('y')
zlabel('z')

