load fisheriris.mat
%get all the classes
class = unique(species);

X = meas;

Y = [];
for ii =1:length(class)
    Y(strcmp([species], class(ii)),:) = ii;
end

%plot em
figure
gscatter3(X(:,1),X(:,2),X(:,3),class(Y)')
title('data with true class')

%%

% format output
[YBinary,cslist,mapping] = NNFormatOutput(Y);

%intialize
layers = 2;
nodes = [4,3];
[NN] = NNInitialize(X,YBinary,layers,nodes);


alpha = 0.6; %learning rate (0 - 1)
max_epoch = 10000;
accuracy_target = 0.05;
%train
[NN,MSE] = NNTrain(X,YBinary,NN,alpha,max_epoch,accuracy_target);

%% Check classifiying

%online NN classify
Yestb = NNOnline(NN, X);

%convert back to our class labels
Yest = NNUnformatOutput(Yestb,cslist,mapping);
    
figure
gscatter3(X(:,1),X(:,2),X(:,3),class(Yest)')
title('data with est class')