% Example (simple Gaussian mixture)
clear all

fontS = 13;

nn = 500;
set1 = randn(nn,2); %around 0,0
set2 = [randn(nn,1).*2+11, randn(nn,1).*2]; %around 5,0
set3 = [randn(nn,1).*2+5, randn(nn,1).*2]; %around 3,-3

% generate a Gaussian mixture with distant modes
data=[set1;set2;set3]; %place all sets in (nn*numsets) x 2 vector 

% call the routine
[bandwidth1,density1,X1,Y1]=kde2d(set2);
[bandwidth2,density2,X2,Y2]=kde2d(set3);
% plot the data and the density estimate
fprintf('bandwidth is %f \n',bandwidth1);

figure(1);
contour3(X1,Y1,density1,50);
hold on
contour3(X2,Y2,density2,50);
colormap Winter
hold on
h1 = plot(set2(:,1),set2(:,2),'r.','MarkerSize',5);
h2 = plot(set3(:,1),set3(:,2),'g.','MarkerSize',5);
hold off

title('Probability Distributions for Two Classes','FontSize',fontS)
xlabel('X1','FontSize',fontS)
ylabel('X2','FontSize',fontS)
zlabel('Probability','FontSize',fontS)
legend([h1,h2],'Class 1 Data','Class 2 Data')
axis([-3 18 -6 6 0 0.04])
