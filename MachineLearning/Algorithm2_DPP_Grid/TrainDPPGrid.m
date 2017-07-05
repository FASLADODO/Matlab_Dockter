function [Model] = TrainDPPGrid(Data,Labels,split)

    %number of classes and data sizes
    cslist = unique(Labels);
    [NN,SS] = size(Data);


    % find limits
    mins = min(Data,[],1);
    maxs = max(Data,[],1);

    % get grid
    ticks = round( (maxs - mins)/split);
    grid = [];
    for ii = 1:SS
        grid(:,ii) = linspace(mins(ii),maxs(ii),split+1);
    end
    
    %get indices (so that one for loop can index into n-dimensional grid)
    boundsidx = [ones(1,SS);ones(1,SS)*size(grid,1)-1];
    gridindices = ndimgrid(boundsidx,split);

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
    Model.dir = [-1,1];
    Model.indices = gridindices;
    Model.f_classify = @sign; %@round
    
    %allocate our space
    ProbabilityRegion = [];
    DensityRegion = [];
    KLRegion = [];
    DPP_Grid_KL = zeros(ones(1,SS)*split); %only for online data
    DPP_Grid_OK = zeros(ones(1,SS)*split); %only for online data
    DPP_Grid_Density = zeros([ones(1,SS)*split]); %only for online data
    DPP_Grid_Probability = zeros([ones(1,SS)*split,length(cslist)]); %only for online data
    DPP_Grid_SWeights = zeros(ones(1,SS)*split);
    GridSize = size(DPP_Grid_KL);
    
    % find data in each region
    for ii = 1:size(Model.limits,1)
        %get limits
        minz = reshape(limits(ii,1,:),1,SS);
        maxz = reshape(limits(ii,2,:),1,SS);

        %segment region
        mask=bsxfun(@ge,Data,minz)& bsxfun(@le,Data,maxz); %check if all dimensions are greater than min and less than max
        i1 = all(mask,2);
        regiondata = Data(i1,:);
        regionlabels = Labels(i1,:);

        %how many points are there
        density = size(regiondata,1);
        probs = [];
        nclass = zeros(2,1);
        
        if(density == 0 )
            %dont store if no data
            probs = ones(length(cslist),1)/length(cslist);
        else
            %find which classes are in region
            for cc = 1:length(cslist)
              %find all indices for the current class
              classidx = regionlabels(:,1) == cslist(cc);
              %grab data for this class in this region
              classregion = regiondata(classidx,:);
              nclass(cc) = size(classregion,1) + 1; %min nclass = 1
            end
            
            %compute probability it this region
            probs = nclass ./ density;
        end
        
        %compute KL divergence
        %THESE ARE THE SAME THING BECAUSE SUM(PROBS)=1 DUMMY
        KL = log10(probs(2)/ probs(1));
        %KL = (probs(2)+probs(1)).*log10(probs(2)./probs(1));
        

        %store in array
        ProbabilityRegion(ii,:) = probs;
        DensityRegion(ii,:) = density;
        KLRegion(ii,:) = KL;
        
        %Putting the grid together
        %figure out which grid we are in
        dist2org = Model.means(ii,:) - Model.origin;
        indexgrid = floor(dist2org ./ Model.steps) + 1;
        selectElement = num2cell(indexgrid);
        indx=sub2ind(GridSize,selectElement{:});
        
        %get current index the easy way
%         indexgrid = Model.indices(ii,:);
%         selectElement = num2cell(indexgrid);
%         indx=sub2ind(GridSize,selectElement{:});
        
        %store in mutlidimensional grid
        DPP_Grid_KL(indx) = KL;
        DPP_Grid_Probability(selectElement{:},:) = probs;
        DPP_Grid_Density(indx) = density;
        DPP_Grid_OK(indx) = density > 1;
    end
    %scale the density 0-1
    DensityRegion = DensityRegion ./ max(DensityRegion(:));
    DPP_Grid_Density = DPP_Grid_Density ./ max(DPP_Grid_Density(:));
    
    %Compute S_{j,k} weights
    DPP_Grid_SWeights = DPP_Grid_KL .* DPP_Grid_Density;
    
    %store
    Model.ProbabilityRegion = ProbabilityRegion;
    Model.DensityRegion = DensityRegion;
    Model.KLRegion = KLRegion;
    Model.SWeightRegion = KLRegion.*DensityRegion;
    Model.DPP_Grid_KL = DPP_Grid_KL;
    Model.DPP_Grid_Probability = DPP_Grid_Probability;
    Model.DPP_Grid_Density = DPP_Grid_Density;
    Model.DPP_Grid_OK = DPP_Grid_OK;
    Model.DPP_Grid_SWeights = DPP_Grid_SWeights;
end




