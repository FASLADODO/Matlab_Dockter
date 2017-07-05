nn = 50;

x1 = randn(nn,1)+3;
x2 = randn(nn,1)*0.6;

X = [x1;x2];
Labels = [ones(nn,1)*1;ones(nn,1)*2];


plotstr = {'Novice';'Expert'};
fsize = 14;
figure
boxplot(X,plotstr(Labels))
xlabel('Skill Level','FontSize',fsize);
ylabel('Segment Counts','FontSize',fsize);
title('Segment Counts Outsize True Expert Region','FontSize',fsize);