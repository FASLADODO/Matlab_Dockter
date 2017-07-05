%timing results for algorithms

numdatapointstrain = 200200;
numdatapointsonline = 1000;
numdimensions = 3;

%dls
dlsiters = 100;
DLStraintime = [9.990358, 9.967117, 9.949694];
DLSonlinetime = [0.008516, 0.008265, 0.008561];

%dpp
dppiters = 10;
DPPtraintime = [189.171356, 198.817361, 189.616518];
DPPonlinetime = [0.479807, 0.469848, 0.467125 ];

%reliefrbf
rbfiters = 10;
rbftraintime = [115.887458, 114.235953, 115.230325];
rbfonlinetime = [67.322508,66.079682,65.864348];

%random forests
rfiters = 10;
rftraintime = [3647.738358,3600.444649, 3829.246879];
rfonlinetime = [15.904296,14.934093,15.064090];

%random forests
nniters = 10;
nntraintime = [87.711817,60.290500,250.788060];
nnonlinetime = [0.100672,0.095283,0.094778];

%% get mean times

dlsmeantrain = mean(DLStraintime) / dlsiters;
dppmeantrain = mean(DPPtraintime) / dppiters;
rbfmeantrain = mean(rbftraintime) / rbfiters;
rfmeantrain = mean(rftraintime) / rfiters;
nnmeantrain = mean(nntraintime) / nniters;


dlsmeanonline = mean(DLSonlinetime) / dlsiters;
dppmeanonline = mean(DPPonlinetime) / dppiters;
rbfmeanonline = mean(rbfonlinetime) / rbfiters;
rfmeanonline = mean(rfonlinetime) / rfiters;
nnmeanonline = mean(nnonlinetime) / nniters;


%% training time
fsize = 14;

c = {'DLS','DPP','RELIEF-RBF','Random Forest','Neural Network'};
y = [dlsmeantrain,dppmeantrain,rbfmeantrain,rfmeantrain,nnmeantrain];
figure
barh(y)
set(gca,'XScale','log')
set(gca,'FontSize',fsize,'yticklabel',c);
xlabel('Time (s)','FontSize',fsize);
ylabel('Algorithm','FontSize',fsize);

%% online time
fsize = 14;

c = {'DLS','DPP','RELIEF-RBF','Random Forest','Neural Network'};
y = [dlsmeanonline,dppmeanonline,rbfmeanonline,rfmeanonline,nnmeanonline];
figure
barh(y)
set(gca,'XScale','log')
set(gca,'FontSize',fsize,'yticklabel',c);
xlabel('Time (s)','FontSize',fsize);
ylabel('Algorithm','FontSize',fsize);



