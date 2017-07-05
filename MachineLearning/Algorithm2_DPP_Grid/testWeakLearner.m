%first some super simple data

Data = [1,1;
        3,4;
        3,6;
        6,3;
        8,5;
        5,5;
        5,6];
Labels = [1;1;1;1;2;2;2];

figure
gscatter(Data(:,1),Data(:,2),Labels)

[NN,SS] = size(Data);
cslist = unique(Labels);

dd = 2;
dtemp = Data(:,dd);
W = ones(NN,1) / NN;
[Thresh,Err,Dir,H_x] = WeakLearner(dtemp,W,Labels,cslist);

Err
Thresh
H_x

%% now some fancier shit
%test weak learner

nn = 20;
scale = 0.3;

d1 = randn(nn,2)*scale;
d2 = randn(nn,2)*scale + repmat([0,1],nn,1);
d3 = randn(nn,2)*scale + repmat([1,0],nn,1);
d4 = randn(nn,2)*scale + repmat([1,1],nn,1);

Data = [d1;d2;d3;d4];
class1 = 2;
class2 = 1;
Labels = [ones(nn,1)*class1; ones(nn,1)*class2; ones(nn,1)*class2; ones(nn,1)*class2];

figure
gscatter(Data(:,1),Data(:,2),Labels)

%% try it on weak

[NN,SS] = size(Data);
cslist = unique(Labels);

errs = [];
thresholds = [];
dirs = [];
%initialize weights
W = ones(NN,1) / NN;
for dd = 1:SS
    dtemp = Data(:,dd);
    [Thresh,Acc,Dir] = WeakLearner(dtemp,Labels,cslist);
    thresholds{dd} = Thresh;
    accs{dd} = Acc;
    dirs{dd} = Dir;
    
end

Bounds = DataBounds(Data);

figure 
gscatter(Data(:,1),Data(:,2),Labels)
hold on
dd = 1;
plot([thresholds{dd};thresholds{dd}],[Bounds(1,dd);Bounds(2,dd)],'r')
hold on
dd = 2;
plot([Bounds(1,dd);Bounds(2,dd)],[thresholds{dd};thresholds{dd}],'b')
hold off

%% Try classifying

storeWeight = [];
for tt = 1:NN
    dtemp = Data(tt,:);
    for dd = 1:SS
        isbelow = dtemp(:,dd) < thresholds{dd};
        classest(dd) = dirs{dd}(isbelow + 1);
        weights(dd) = accs{dd};
    end
    storeWeight(tt,:) = sum(classest.*weights) / sum(weights);
    LabelsEst(tt,:) = round(sum(classest.*weights) / sum(weights) );
    
end

errall = mean(LabelsEst ~= Labels)
    
figure 
gscatter(Data(:,1),Data(:,2),LabelsEst)
    
figure 
scatter(Data(:,1),Data(:,2),10,storeWeight)
colorbar
colormap cool
    











