% test the subtask partitioning (like nisky group)

%% first create rand data

mu1 = [3,3,3.6];
sc1 = [2,0.1,0.01];

mu2 = [-5,8,4];
sc2 = [0.01,5,0.01];

SS = length(sc1);
nn = 200;

d1 = randn(nn,SS).*repmat(sc1,nn,1) + repmat(mu1,nn,1);
d2 = randn(nn,SS).*repmat(sc2,nn,1) + repmat(mu2,nn,1);

D = [d1;d2];
labels = [ones(nn,1)*-1; ones(nn,1)*1];

figure
gscatter3(D(:,1),D(:,2),D(:,3),labels)
xlabel('x')
ylabel('y')
zlabel('z')

%% Now try kmeans to label each sub motion (ignore the states)


