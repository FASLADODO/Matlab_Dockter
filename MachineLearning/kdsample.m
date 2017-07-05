nn = 100;

Data = randn(nn,2)*3;



kn = 10;

Mdl{1} = KDTreeSearcher(Data);
Mdl{2} = KDTreeSearcher(Data);

testpoint = Data(25,:);
idwd = knnsearch(Mdl{1},testpoint,'K',kn);
windowData = Data(idwd,:);

figure
scatter(testpoint(:,1),testpoint(:,2),'g*')
hold on
scatter(Data(:,1),Data(:,2),'r.')
hold on
scatter(windowData(:,1),windowData(:,2),'b.')
hold off