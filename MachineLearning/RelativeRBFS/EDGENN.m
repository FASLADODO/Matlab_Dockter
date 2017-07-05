%% SPECIAL VERSION THAT CAN BE RUN IN CITRIX

%load DataTrain, LabelsTrain, DataTest, LabelsTest, DataValidate,
%LabelsValidate, and key

tic
fprintf('loading the saved data...')

load('DataMatrixEdge.mat')

fprintf('done \n')
toc

%% now that we have the data, choose which columns and train the nueral net


usecolumns = [key.col.dx,key.col.dy,key.col.dz,key.col.ddx,key.col.ddy,key.col.ddz,key.col.velalpha,key.col.velbeta,key.col.velmag,key.col.accmag,key.col.accalpha,key.col.accbeta];


Xtrain = DataTrain(:,usecolumns)';
%get NN labels
targets = getNNLabels(LabelsTrain)';

%%

%makes the nnet object
net = patternnet(25);
%view(net)

%train and view
net = train(net,Xtrain,targets);

%% Test net on validation set

Xvalidate = DataValidate(:,usecolumns)';
targetsvalidate  = getNNLabels(LabelsValidate)';

%estimate classes
netout = net(Xvalidate);
perf = perform(net,targetsvalidate,netout);

%check accuracy
[~,classest] = max(netout);
corr = classest == LabelsValidate';
acc = mean(corr)

figure
plot(netout(1,:),'r')
hold on
plot(netout(2,:),'b')
hold off
legend('class 1','class 2')