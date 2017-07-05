function Model = LogRegTrain(Data,Labels)
%LogRegTrain: train logistic regression parameters
% Data = ND matrix with each row as a sample
% Labels = rowwise labels

[NN,SS] = size(Data);
[NL,SL] = size(Labels);
classes = unique(Labels);

if(SL > 1)
   error('give single unique value for each class') 
end


multiclass = 0;
if(length(classes) > 2)
   warning('only two classes supported!') 
   multiclass = 1;
end

if(multiclass)
    %Now set labels to 0 or 1
    P = zeros(NN,length(classes));
    for cc = 1:length(classes)
        temp = zeros(1,length(classes));
        temp(cc) = 1;
        Model.Ref(cc,:) = temp;
        idcs = Labels == classes(cc);
        P(Labels == classes(cc),: ) = repmat(temp,sum(idcs),1);
    end
else
    %Now set labels to 0 or 1
    P = zeros(NN,1);
    for cc = 1:length(classes)
        Model.Ref(cc) = cc-1;
        P(Labels == classes(cc) ) = cc - 1;
    end
end
%get probabilities (scale everything by 0.001;
P = abs(P - eps);

%Get logit
logodds = logit(P);

%Now Try and get params
if(multiclass)
    for ii = 1:length(classes)
        tp = pinv(Data)*logodds(:,ii);
        Model.Params(:,ii) = tp;
    end
else
    Model.Params = pinv(Data)*logodds;
end

end