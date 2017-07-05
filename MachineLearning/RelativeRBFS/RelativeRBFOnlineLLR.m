function [LLR,Class,ProbClass] = RelativeRBFOnlineLLR(Data,Model,Thresh)
%Model comes from RelativeRBFTrain
%Data: samples as rows, dimensions as columns

    if(nargin < 3)
       Thresh = Model{end}.ThresholdSep; %0.05; 
    end
    Prob_All = [];
    columnID = [];
    for cc = 1:length(Model)
        for idc = 1:Model{cc}.TotalClusters;
            ptemp = gaussianEval(Data,Model{cc}.cluster{idc});
            ProbClass{cc}(:,idc) = ptemp;
        end
        Prob_All(:,cc) = max(ProbClass{cc},[],2);
    end

    Prob_All(max(Prob_All,[],2) < Thresh,:)=[];
    
    %[LLR,combos] = cumSumLLR(Prob_All);
    LLR = log(sum(Prob_All(:,1)) / sum(Prob_All(:,2)) );
    
    %take class from column
    Class = (LLR < 0) + 1;
    
end