function [ testClass, arrClass, listClass ] = DPP_Discriminant_Online( tData, oData, region, weights, x1grid, x2grid, subsample )
%Discriminant Phase Portrait Online Classification
% tData is the training data for all classes and all runs
% oData is the online sample data
% kk = number of classes
% ii = number of states
% nn = number of data points in all training data
% tData{kk}.state{ii}() = [1 x nn]
% oData.state{ii}
%region has state data embedded in struct for each region
% weights are the seperability measures
% x1grid is the regions in state 1
% x2grid is the regions in state 2
% subsample is the subsamples rate ie 1 in 10 samples


%NN size (This should be odd for science)
K = 5;

%Prune data (could be done in seperate function)
train = [];
grouper = [];

cc = 1;
for kk = 1:length(tData)
    idx = 1:subsample:length(tData{kk}.state{1});
    train = [train; tData{kk}.state{1}(idx), tData{kk}.state{2}(idx) ];
    grouper = [grouper; ones(length(idx),1)*kk ];
    length(idx)
end

%online data
new = [];
new = [oData.state{1},oData.state{2}];

% sums/arrays for tracking estimate
testClass = zeros(1,length(region));
listClass = [];
arrClass = [];
for tt = 1 : length(new)
    sample = new(tt,:);
    idx1 = find(x1grid > sample(1),1) - 1;
    idx2 = find(x2grid > sample(2),1) - 1;
    if(idx1 < 1)
       idx1 = 1; 
    end
    if(idx2 < 1)
       idx2 = 1; 
    end
    if(idx1 > length(x1grid))
       idx1 = length(x1grid); 
    end
    if(idx2 > length(x2grid))
       idx2 = length(x2grid); 
    end
    
    region_train = [];
    region_group = [];
    haveData = 0;
    for kk = 1:length(region)
        if(~isempty(region{kk}.hor{idx1}.vert{idx2}.state{1}) )
            region_train = [region_train; region{kk}.hor{idx1}.vert{idx2}.state{1}(1:subsample,end), region{kk}.hor{idx1}.vert{idx2}.state{2}(1:subsample,end) ];
            region_group = [region_group; ones(length(region{kk}.hor{idx1}.vert{idx2}.state{1}(1:subsample,end)),1)*kk ];
            haveData = 1;
        end
    end
    
    if(haveData)
        [ IDX,Class,ClassArr,Dist] = KNN_RLD( region_train, region_group, sample, K );

        arrClass(tt,:) = ClassArr;
        for kk = 1:length(region)
           if kk == Class
               listClass(tt) = kk;
               testClass(kk) = testClass(kk) + weights(idx1,idx2); %*(1/(mean(Dist) + 1));
               break
           end
        end
    else
        arrClass(tt,:) = [];
        listClass(tt) = 0; %0 is i dont know
    end
    
    
end



end

