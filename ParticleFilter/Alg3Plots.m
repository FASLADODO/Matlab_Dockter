% Tim function plots

clear all

fontS = 14;

%runs
run = 10;
noise = 0.2;

%params
% Phi1 = [5, -250];
% Phi2 = [-5, 250];
Phi1 = [4, 0];
Phi2 = [1, 0];

%inputs
% x1 = 50:0.1:75;
% x2 = 25:0.1:50;

x1 = 0:0.1:5;
x2 = 0:0.1:5;


hFig = figure(1)
for kk = 1:run
    %noise input
    x1_noise = x1 + (rand(1,length(x1))-.5)*noise;
    x2_noise = x2 + (rand(1,length(x2))-.5)*noise;

    %functions
%     F1_temp = Phi1(1)*x1_noise + Phi1(2);
%     F2_temp = Phi2(1)*x2_noise + Phi2(2);
    F1_temp = Phi1(1)*x1_noise + Phi1(2);
    F2_temp = Phi2(1)*x2_noise.^2 + Phi2(2);
    
    %noise ouput
    F1_temp = F1_temp + (rand(1,length(F1_temp))-.5)*noise;
    F2_temp = F2_temp + (rand(1,length(F2_temp))-.5)*noise;
    
    %store
    F1_s(kk,:) = F1_temp;
    F2_s(kk,:) = F2_temp;
    X1_s(kk,:) = x1_noise;
    X2_s(kk,:) = x2_noise;
    
    %limits
%     idx1 = find(F1_temp<0);
%     idx2 = find(F2_temp<0);
%     F1_temp(idx1)=[]; %0
%     F2_temp(idx2)=[]; %0
%     x1_noise(idx1)=[]; %50 + (rand(1)-.5)*5;
%     x2_noise(idx2)=[]; %50 + (rand(1)-.5)*5;
    
    %plot
    h1 = scatter(x1_noise,F1_temp,'ro');
    hold on 
    h2 = scatter(x2_noise,F2_temp,'gx');
    hold on
    
end
hold on

%reshape arrays
data1 = [ reshape(F1_s,[],1), reshape(X1_s,[],1) ];
data2 = [ reshape(F2_s,[],1), reshape(X2_s,[],1) ];

%get KDE
[bandwidth_1,density_1,X_1,Y_1]=kde2d(data1);
[bandwidth_2,density_2,X_2,Y_2]=kde2d(data2);

contour3(Y_1,X_1,density_1,50);
hold on
contour3(Y_2,X_2,density_2,50);
hold off


grid on
view(0, 45);
xlabel('x','FontSize',fontS)
ylabel('xdot','FontSize',fontS)
zlabel('Probability','FontSize',fontS)
title('State space with probabilities','FontSize',fontS)
legend([h1,h2],'class 1', 'class 2','FontSize',fontS)

axis([0,6,0,30])
