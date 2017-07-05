%test RBF Thesis

rng('default')

%populate some data
fsize = 14

nn = 2000;

mu1 = 2;
sg1 = 0.2;
x1 = [randn(nn,1)*sg1 + mu1, randn(nn,1)*sg1 + mu1 ];

mu2 = 2.3;
sg2 = 0.2;
x2 = [randn(nn,1)*sg2 + mu2, randn(nn,1)*sg2 + mu2 ];

%combine it
Data = [x1; x2];
Labels = [ones(nn,1)*-1; ones(nn,1)*1];

figure
gscatter(Data(:,1),Data(:,2),Labels)

xlabel('x1','FontSize',fsize)
ylabel('x2','FontSize',fsize)
% zlabel('W_{KL}','FontSize',fsize)
[lh,ic,ip,it]=legend('show');
lh.FontSize = fsize;
lh.Location = 'NorthEast';
% title('DPP KL Divergence')
% colormap(flipud(cool))
% hc = colorbar;
% ylabel(hc, 'W_{KL}','FontSize',fsize)

%% Get single class rbfs and plot

[Difference,ClassData,Model] = SimpleRelativeRBFTrain(Data,Labels);

%%Combine it
AllData = [ClassData{1}; ClassData{2}];
AllDiff = [Difference{1}; Difference{2}];

figure
gscatter(Data(:,1),Data(:,2),Labels)
hold on
Surface3D(AllData(:,1),AllData(:,2),AllDiff,'mesh');
hold off
xlabel('x1','FontSize',fsize)
ylabel('x2','FontSize',fsize)
zlabel('W_{RBF}','FontSize',fsize)
[lh,ic,ip,it]=legend('show');
lh.FontSize = fsize;
lh.Location = 'NorthEast';
% title('DPP KL Divergence')
colormap(flipud(cool))
hc = colorbar;
ylabel(hc, 'W_{RBF}','FontSize',fsize)

%% Now we sort our diffs

[Diffz,Dataz,Labelz] = GPRSubsample(Data,Labels);

ratio = 0.5;
classes = length(unique(Labelz));

NN = length(Diffz);

[sortDiff,idd] = sort(Diffz,'descend');

thresh = sortDiff( floor(NN*ratio), :);

figure
p1=plot(sortDiff);
hold on
p2=plot( [1,NN] , [thresh, thresh], 'k');
hold off
xlabel('sample #','FontSize',fsize)
ylabel('W_{i,rbf}','FontSize',fsize)
lh=legend('W_{T}','T_{sub}');
lh.FontSize = fsize;
p1(1).LineWidth = 2;
p2(1).LineWidth = 2;

%% subsample that data

%n*r_{w}
lim = floor(NN*ratio);

%subsamples
subdiff = sortDiff(1:lim);
subdata = Dataz(idd(1:lim),:);
sublabels = Labelz(idd(1:lim),:);

figure
p=gscatter(subdata(:,1),subdata(:,2),sublabels);

xlabel('x1','FontSize',fsize)
ylabel('x2','FontSize',fsize)
% zlabel('W_{KL}','FontSize',fsize)
[lh,ic,ip,it]=legend('show');
lh.FontSize = fsize;
lh.Location = 'NorthEast';

%% lets train gpr with the subsample data

GPRMDL = GPRClassifyTrain(subdata,sublabels);

[Yest,sigmaest] = GPRClassifyOnline(GPRMDL,Data);

ClassEst = sign(Yest);

figure
scatter(Data(:,1),Data(:,2),10,Yest);
xlabel('x1','FontSize',fsize)
ylabel('x2','FontSize',fsize)
colormap(flipud(cool))
hc = colorbar;
ylabel(hc, '\mu_{*}','FontSize',fsize)

figure
gscatter(Data(:,1),Data(:,2),ClassEst);
xlabel('x1','FontSize',fsize)
ylabel('x2','FontSize',fsize)
[lh,ic,ip,it]=legend('show');
lh.FontSize = fsize;
lh.Location = 'NorthEast';

corr = Labels == ClassEst;
acc = mean(corr)