% Nearest Neighbor

fontS = 14;

dsize1 = 30;
dsize2 = 45;
% train data
data1 = [randn(1,dsize1) + 3; randn(1,dsize1) + 3];
data2 = [randn(1,dsize2) + 1; randn(1,dsize2) + 1];


%new point
nsize = 3;
new = [randn(1,nsize) + 2; randn(1,nsize) + 2];

figure(1)
plot(data1(1,:),data1(2,:),'go','LineWidth',2,'MarkerSize',10);
hold on
plot(data2(1,:),data2(2,:),'bx','LineWidth',2,'MarkerSize',10);
hold on
plot(new(1,:),new(2,:),'ro','MarkerFaceColor',[1 0 0],'MarkerSize',10);
hold on

%number of neighbors
K = 4;

%combine data
totaldata = [data1, data2];
%groups
grouper = [ones(dsize1,1);ones(dsize2,1)*2];

%classify and get indexes
class = knnclassify(new', totaldata', grouper')
IDX = knnsearch(totaldata',new','K',K)

%plot indexs
plot(totaldata(1,IDX),totaldata(2,IDX),'kp','MarkerSize',12);

hold off
xlabel('x position','FontSize',fontS)
ylabel('y position','FontSize',fontS)
title('Nearest Neighbor (K = 4)','FontSize',fontS)
legend('Training Class 1','Training Class 2','New data','NN','FontSize',fontS)

%% 

[ IDX,Class,ClassArr,Dist] = KNN_RLD( totaldata', grouper', new(:,1)', K )