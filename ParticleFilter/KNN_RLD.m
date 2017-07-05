function [ IDX,Class,ClassArr,Dist ] = KNN_RLD( train, group, sample, K )
%perform k nearest neighbors search using trainign data
% train = ss x n (training data)
% group = 1 x n (class labels)
% sample = ss x m (new incoming data)
% K = number of nearest neighbors

%Check sizes
[row1,col1] = size(train);
[row2,col2] = size(sample);

if(col1 ~= col2)
   error('sample and training must be same size') 
end

%initialize
Dist = zeros(row2,K);
ClassArr = zeros(row2,K);
Class = zeros(1,row2);
IDX = zeros(row2,K);

% loop through all samples
for kk = 1:row2

    minNN = ones(1,K)*10000;
    minIDX = ones(1,K);
    minClass = ones(1,K);
    for jj = 1:row1
        %compute distance
        dist = 0;
        for ss=1:col2
            dist = dist + (sample(kk,ss) - train(jj,ss) )^2;
        end
        dist = sqrt(dist);
        
        %get indexes and shift
        idx = find(dist<minNN,1);
        if(~isempty(idx))
            for ii = K:-1:idx+1
                minNN(ii) = minNN(ii-1);
                minClass(ii) = minClass(ii-1);
                minIDX(ii) = minIDX(ii-1);
            end
            minNN(idx) = dist;
            minClass(idx) = group(jj);
            minIDX(idx) = jj;
        end
    end

    %save NN for each sample
    Dist(kk,:) = minNN;
    ClassArr(kk,:) = minClass;
    IDX(kk,:) = minIDX;
    Class(kk) = mode(minClass);
end

