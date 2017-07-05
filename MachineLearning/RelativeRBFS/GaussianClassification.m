function [Class,Accuracy,Probability] = GaussianClassification(Data,Labels)
%Data = class data to train on (rows are samples, columns are dimensions)
%Labels = column of unique class lables with same length as Data [1;2...]

    [NN,SS] = size(Data);
    cslist = unique(Labels);

    Prob_All = [];
    for cc = 1:length(cslist)
        Dtemp = Data(Labels == cslist(cc),:);
        mu = mean(Dtemp);
        sigma = cov(Dtemp);
        ptemp = mvnpdf(Data,mu,sigma);

        Prob_All(:,cc) = ptemp;
    end

    %get the highest probability for each data points
    [Probability,maxid] = max(Prob_All,[],2);
    
    %take class from column
    Class = cslist(maxid);

    %get classifcation accuracy
    Correct = Class == Labels;

    %get accuracy
    Accuracy = mean(Correct);

end