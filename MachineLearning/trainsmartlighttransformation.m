% find new transformation for smart light pixel data

DNew = load('C:/temp/pixelpointsnew.csv');

DOld = load('C:/temp/camerapoints.csv');


rows = ones(length(DOld),1)*640;
camframe = DOld;
camframe(:,2) = rows - camframe(:,2);

figure
scatter(camframe(:,1),camframe(:,2),500,'b.')
axis([0 640 0 480])

figure
scatter(DNew(:,1),DNew(:,2),500,'r.')
hold on
scatter(DOld(:,1),DOld(:,2),500,'b.')
hold off
axis([0 640 0 480])

rnmin = min(DNew);
rnmax = max(DNew);
romin = min(DOld);
romax = max(DOld);

sn = rnmax - rnmin
so = romax - romin

scale = mean(sn ./ so)

scaleinc = 1/scale

%% old to new

TTest = Transformation2D(0.43,[-304,83],0.629)
% TTest = Transformation2D(0,[0,0],0.65)

DOld_temp = [DOld, ones(length(DOld),1)];

DOld_tran = TTest * DOld_temp' ;

DOld_tran = DOld_tran';

DOld_tran = DOld_tran./repmat(DOld_tran(:,3),1,3);



figure(1)
scatter(DNew(:,1),DNew(:,2),1000,'r.')
hold on
scatter(DOld_tran(:,1),DOld_tran(:,2),1000,'b.')
hold off
legend('NewPoints','OldPoints Transform','location','southwest')
axis square

%% save it

dlmwrite('T_Prior.rod',TTest)