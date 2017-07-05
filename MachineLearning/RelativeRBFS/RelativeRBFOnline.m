function [Class,Probability] = RelativeRBFOnline(Data,Model,Thresh)
%Model comes from RelativeRBFTrain
%Data: samples as rows, dimensions as columns

    if(nargin < 3)
       Thresh = Model{end}.ThresholdSep; %0.05; 
    end
    Prob_All = [];
    columnID = [];
    idx = 1;
    for cc = 1:length(Model)
        for idc = 1:Model{cc}.TotalClusters;
            ptemp = gaussianEval(Data,Model{cc}.cluster{idc});
            Prob_All(:,idx) = ptemp;
            columnID(idx) = cc;
            idx = idx + 1;
        end
    end
    %get the highest probability for each data points
    [Probability,maxid] = max(Prob_All,[],2);
    
    %take class from column
    Class = columnID(maxid)';
    
    %Only keep class estimate if probability is high enough
    Class(Probability < Thresh) = 0;
    
end