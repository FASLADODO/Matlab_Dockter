function [idxtrain, idxval, idxtest] = dividerand(NN,trainRatio,valRatio,testRatio) 
%should be a functional match to the one in NN net lib
%Get Training, Validation, and Test data from main data set
%trainRatio = 0.7, valRatio = 0.15, testRatio = 0.15

%[NN,~] = size(Data);
%check sums
if(abs(sum([trainRatio,valRatio,testRatio]) - 1) > eps)
   error('ratios must sum to one');
end

%get number of points to sample
kktr = round(trainRatio*NN);
kkva = round(valRatio*NN);
kkte = round(testRatio*NN);

%make sure we didn't mess up our fractions
while(kktr + kkva + kkte < NN)
   kktr = kktr + 1; 
end
while(kktr + kkva + kkte > NN)
   kktr = kktr - 1; 
end

%get ourstarting vector
rng('shuffle')
AllSamples = randperm(NN);

%get sub sample indices from all samples (this prevents replacement)
idxtrain = sort( AllSamples(1:kktr) );
idxval = sort( AllSamples(kktr+1:kktr+kkva) );
idxtest = sort( AllSamples(kktr+kkva+1:kktr+kkva+kkte) );


end