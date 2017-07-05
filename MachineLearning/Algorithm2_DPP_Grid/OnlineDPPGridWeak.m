function [ClassStore,RawStore] = OnlineDPPGridWeak(Data,Model)
%Do online classification for time series Data using Model
%Data: rows as samples, #columns should match that trained with Model
%Model: DPP grid model trained in TrainDPPGrid.m

[NN,SS] = size(Data);

ClassStore = [];
RawStore = [];
%go through each timestamp
for tt= 1:NN

    %see if we're in bounds
    if(any(bsxfun(@lt,Data(tt,:),Model.origin),2) || any(bsxfun(@ge,Data(tt,:),Model.maxlimit),2) )
        ClassAll = 0;
        RawSum = 0;
    else
        %if so figure out which grid region we're in
        dist2org = Data(tt,:) - Model.origin;
        indexgrid = floor(dist2org ./ Model.steps) + 1;
        selectElement = num2cell(indexgrid);
        indx=sub2ind(size(Model.DPP_Grid_Weight),selectElement{:});
        %get weights and thresholds for that region
        Weight_temp = Model.DPP_Grid_Weight_All(selectElement{:},:);
        Thresh_temp = Model.DPP_Grid_Thresh(selectElement{:},:);
        Class_temp = Model.DPP_Grid_Class(selectElement{:},:,:);
        UsePT = Model.DPP_Grid_OK(selectElement{:},:);
        
        %make sure we know something about this region
        if(UsePT ~= 1)
            ClassAll = 0;
            RawSum = 0;
        else
            %we're good, classify
            %loop through each dimension
            classest = [];
            weights = [];
            for dd = 1:SS
                isbelow = Data(tt,dd) < Thresh_temp(dd);
                dir = Class_temp(:,:,dd,:);
                classest(dd) = dir(isbelow + 1);
                weights(dd) = Weight_temp(dd);
            end

            %Store Classification
            RawSum = sum(classest.*weights) / sum(weights);
            ClassAll = Model.f_classify( RawSum );
        end
    end
    
    %store it for this point
    ClassStore(tt,:) = ClassAll;
    RawStore(tt,:) = RawSum;
end


end