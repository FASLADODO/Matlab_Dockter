function [Model] = GaussianClassifierTrain(Data,Labels)
%Taken from CSCI 5525 Leg04-gendisc

    [NN,SS] = size(Data);
    gc = unique(Labels);

    SIGMA = zeros(SS,SS);
    for cc = 1:length(gc)
        temp = Data(Labels == gc(cc),:);
        PC{cc} = length(temp)/NN; %class conditionals

        GM{cc}.Sigma = cov(temp);
        GM{cc}.Mu = mean(temp);
        
        SIGMA = SIGMA + GM{cc}.Sigma;
    end
    SIGMA = SIGMA ./ length(gc);

    Model.W = inv(SIGMA)*(GM{1}.Mu - GM{2}.Mu)';
    Model.w0 = (-1/2)*GM{1}.Mu*inv(SIGMA)*GM{1}.Mu' + (1/2)*GM{2}.Mu*inv(SIGMA)*GM{2}.Mu' + log10(PC{1}/PC{2});

end