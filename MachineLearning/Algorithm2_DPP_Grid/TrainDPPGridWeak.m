function [Model] = TrainDPPGridWeak(Data,Labels,split)
%train each grid with a weak learner in each dimesnions
%right now weak learner is a decision stump

    %number of classes and data sizes
    cslist = unique(Labels);
    [NN,SS] = size(Data);
    
    %Convert labels to +/- 1 (too make life easier)

    % find limits
    mins = min(Data,[],1);
    maxs = max(Data,[],1);

    % get grid
    ticks = (maxs - mins)/split;
    grid = [];
    for ii = 1:SS
        grid(:,ii) = linspace(mins(ii),maxs(ii),split+1);
    end

    %get grid origin, stepsize and limits
    origin = grid(1,:);
    maxlimit = grid(end,:);
    steps = mean(diff(grid));
    limits = SegmentLimitGrid(grid); %this is a fancy function to get ndimensional grid limtis

    clear MeansAll
    %get means
    for ii = 1:size(limits,1)
        minz = reshape(limits(ii,1,:),1,SS);
        maxz = reshape(limits(ii,2,:),1,SS);
        middle = (maxz + minz) / 2;
        MeansAll(ii,:) = middle;
    end

    % store
    Model.grid = grid;
    Model.origin = origin;
    Model.maxlimit = maxlimit;
    Model.steps = steps;
    Model.limits = limits;
    Model.means = MeansAll;
    Model.dir = cslist;
    Model.f_classify = @sign; %@round

    % find data in each region
    ThresholdRegion = [];
    ResultantRegion = [];
    WeightRegion = [];
    DPP_Grid = zeros(ones(1,SS)*split); %only for online data
    DPP_Grid_OK = zeros(ones(1,SS)*split); %only for online data
    DPP_Grid_Thresh = zeros([ones(1,SS)*split,SS]); %only for online data
    DPP_Grid_Class = zeros([ones(1,SS)*split,SS,length(cslist)]); %only for online data
    DPP_Grid_Weights = zeros([ones(1,SS)*split,SS]); %only for online data
    for ii = 1:size(Model.limits,1)
        %get limits
        minz = reshape(limits(ii,1,:),1,SS);
        maxz = reshape(limits(ii,2,:),1,SS);

        %segment region
        mask=bsxfun(@ge,Data,minz)& bsxfun(@le,Data,maxz); %check if all dimensions are greater than min and less than max
        i1 = all(mask,2);
        regiondata = Data(i1,:);
        regionlabels = Labels(i1,:);
        RegionOK = 0;
        
        if(isempty(regiondata) )
            %dont store
            Thresholds = zeros(1,SS);
            ClassResult = zeros(2,SS); %if below thresh_temp
            W_Store = zeros(1,SS);
        else
            RegionOK = 1;
            Thresholds = [];
            ClassResult = [];
            W_Store = [];
            %for each dimension, pick a weak learner
            for dd = 1:SS
                dtemp = regiondata(:,dd);
                %train weak learner
                [thresh_temp,ACC_temp,Dir_temp] = WeakLearner(dtemp,regionlabels,cslist);

                %store it
                Thresholds(dd) = thresh_temp;
                ClassResult(dd,:) = Dir_temp; %if below thresh_temp
                W_Store(dd) = ACC_temp;
            end
        end
        
        %Putting the grid together
        %figure out which grid we are in
        dist2org = Model.means(ii,:) - Model.origin;
        indexgrid = floor(dist2org ./ Model.steps) + 1;
        selectElement = num2cell(indexgrid);
        indx=sub2ind(size(DPP_Grid),selectElement{:});
        DPP_Grid(indx) = mean(W_Store);
        DPP_Grid_Weights(selectElement{:},:) = W_Store;
        DPP_Grid_Thresh(selectElement{:},:) = Thresholds;
        DPP_Grid_Class(selectElement{:},:,:) = ClassResult;
        DPP_Grid_OK(selectElement{:},:,:) = RegionOK;
        
        %store
        ThresholdRegion(ii,:) = Thresholds;
        ResultantRegion(ii,:,:) = ClassResult;
        WeightRegion(ii,:) = W_Store;
    end

    %store
    Model.ThresholdRegion = ThresholdRegion;
    Model.ResultantRegion = ResultantRegion;
    Model.WeightRegion = WeightRegion;
    Model.DPP_Grid_Weight = DPP_Grid;
    Model.DPP_Grid_Weight_All = DPP_Grid_Weights;
    Model.DPP_Grid_Thresh = DPP_Grid_Thresh;
    Model.DPP_Grid_Class = DPP_Grid_Class;
    Model.DPP_Grid_OK = DPP_Grid_OK;
end