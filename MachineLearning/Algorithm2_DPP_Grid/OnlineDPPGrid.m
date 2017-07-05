function [Class,ClassTime,ScoreTime,ScoreAll] = OnlineDPPGrid(Data,Model)
%Do online classification for time series Data using Model
% Inputs:
% Data: rows as samples, #columns should match that trained with Model
% Model: DPP grid model trained in TrainDPPGrid.m
% Returns:
% Class is the overall class estimate for the segment
% ClassTime is the per timestep class estimate
% ScoreTime is the running sum score at each time step

Score = 0;
ScoreAll = [];
%go through each timestamp
for tt= 1:size(Data,1)

    %see if we're in bounds
    if(any(bsxfun(@lt,Data(tt,:),Model.origin),2) || any(bsxfun(@ge,Data(tt,:),Model.maxlimit),2) )
        pval = 0;
    else
        %if so figure out KL value in that grid region
        dist2org = Data(tt,:) - Model.origin;
        indexgrid = floor(dist2org ./ Model.steps) + 1;
        selectElement = num2cell(indexgrid);
        indx=sub2ind(size(Model.DPP_Grid_SWeights),selectElement{:});
        pval = Model.DPP_Grid_SWeights(indx);
    end
    
    %keep track of overall grasp val
    ScoreAll(tt,:) = pval;
end

%figure out per time class and score
ScoreTime = cumsum(ScoreAll);
ClassTime = Model.f_classify(ScoreTime);
Score = ScoreTime(end);

%figure out the overall class
Class = Model.f_classify(Score);

end