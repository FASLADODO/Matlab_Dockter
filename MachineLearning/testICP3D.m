nn = 100;

x1 = linspace(0,2,nn)';
modelz = x1*2 + 4 + 0.1 + 0.01*randn(nn,1);
modeldata = [x1,x1,modelz];

x2 = linspace(0,2,nn)';
targetz = x2*4 + 4 + 0.1 + 0.01*randn(nn,1);
targetdata = [x2,x2,targetz];


figure
scatter3(modeldata(:,1),modeldata(:,2),modeldata(:,3),'b+')
hold on
scatter3(targetdata(:,1),targetdata(:,2),targetdata(:,3),'ro')
hold off

[Correspondence,Errors,T] = IterativeClosestPoint3D(modeldata,targetdata);

T


figure
scatter3(modeldata(:,1),modeldata(:,2),modeldata(:,3),'b+')
hold on
scatter3(targetdata(:,1),targetdata(:,2),targetdata(:,3),'ro')
hold on
for ii = 1:nn
   plot3([targetdata(ii,1),modeldata(Correspondence(ii),1)],[targetdata(ii,2),modeldata(Correspondence(ii),2)],[targetdata(ii,3),modeldata(Correspondence(ii),3)])
   hold on
end
hold off

transformdata = (T*[modeldata,ones(nn,1)]')';
transformdata = transformdata(:,[1:3]);
figure
scatter3(transformdata(:,1),transformdata(:,2),transformdata(:,3),'b+')
hold on
scatter3(targetdata(:,1),targetdata(:,2),targetdata(:,3),'ro')
hold on
