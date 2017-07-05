function [Model] = RelativeRBFTrain(Data,Labels,kmeanlimit,ploton)
%Data = class data to train on (rows are samples, columns are dimensions)
%Labels = column of unique class lables with same length as Data [1;2...]
%kmeanlimit = estimate of max number of clusters for each class
%ploton = 'ploton' or 'plotoff'

[NN,SS] = size(Data);
fsize = 14;

%optimal bandwidth (rule of thumb)
sig = norm(std(Data));
n = length(Data);
bw = 1.06*sig*(n^(-1/5))

%get all unique classes
if(isnumeric(Labels))
    cslist = unique(Labels);
    cslookup = strread(num2str(cslist'),'%s');
else
    cslookup = unique(Labels);
    cslist = 1:length(cslookup);
end

%threshold
% thresh = log10(2)/(2^SS)
% thresh = (log10(2)/2)*10^(-(SS+1))
% thresh = 2*10^(-((pi/2)*SS));
thresh = 0.1;
% thresh = 2
plotsize3 = 10;

mr = max(cslist) + 1; % should be 3 is cslist = [1,2]
%loop through all classes
Data_All = [];
Diff_All = [];
Class_All = [];
Check_Thresh = [];

%Find all relative RBFS for class data
for cc = 1:length(cslist)
    %current class and opposite class
    DON = Data(Labels == cslist(cc),:);
    DOFF = Data(Labels ~= cslist(cc),:);
    
    nn = size(DON,1);
    
    %compute dem rbfs
    [pwithin,pbetween] = ComputeDiscriminateRBF(DON,DOFF,bw);
   
    Difference{cc} = ComputeRBFDifference(pwithin,pbetween);
    %for seperability plots
    SavePR{cc}{cc} = Difference{cc};
    SavePR{cc}{mr-cc} = ComputeRBFDifference(pbetween,pwithin);
    DataClass{cc} = DON;
    
    %extra stuff
    Data_All = [Data_All; DON];
    Diff_All = [ Diff_All; Difference{cc}];
    Class_All = [Class_All; ones(nn,1)*cc ];
    
    %To make sure our threshold is valid
    Check_Thresh = [ Check_Thresh, max(Difference{cc}) ];
end


%make sure our threshold will work, if not use 1/10th of max
if(max(Check_Thresh) < thresh)
   warning('Data is not very seperable')
   thresh = max(Check_Thresh) / 5
end

%Test Accuracy with Seperability Values
Class_Estimate = [];
for cc = 1:length(cslist)
    %get the current class data and size
    lkp = [cc,mr-cc];
    
%     Dtemp = DataClass{cc};
%     [PROnline] = ComputeDiscriminantRBFOnline(Dtemp,DataClass,bw);
    PROnline = SavePR{cc};
    %both difference probabilities
    COMB = [PROnline{cc},PROnline{mr-cc}];
    
    %find the max ratio
    [DFIN,DID] = max(COMB,[],2);
    %convert to actual class
    ClassG = lkp(DID)';
    %only use classification if above a threshold
    isgood = (DFIN>thresh);
    estc = isgood.*ClassG;
    Class_Estimate = [Class_Estimate; estc];
    ng = size(estc,1);
    %record accuracy
    corr = estc == ones(ng,1)*cc;
    AccuracySeperable(cc) = mean(corr(isgood));

end


%plot diffs
if(strcmp(ploton,'ploton'))
    if(SS == 1)
        figure
        if(length(unique(Class_Estimate)) == 2)
            gclrz = 'br';
        else
            gclrz = 'gbr';
        end
        gscatter(Data_All,zeros(length(Data_All),1),Class_Estimate,gclrz);
        for cc = 1:length(cslist)
            hold on
            scatter(DataClass{cc},Difference{cc},'g.');
        end
        hold off
        str = sprintf('Seperability w/ Discriminant (%f)',mean(AccuracySeperable) );
        title(str)
        xlabel('sample#')
        ylabel('Relative P')
    elseif(SS == 2)
        figure
        if(length(unique(Class_Estimate)) == 2)
            gclrz = 'br';
        else
            gclrz = 'gbr';
        end
        %gscatter(Data_All(:,1),Data_All(:,2),Class_Estimate,gclrz);
        strlab = {'Class 1';'Class 2'};
        gscatter(Data_All(:,1),Data_All(:,2),strlab(Labels),'br');
        xlabel('State 1','FontSize',fsize)
        ylabel('State 2','FontSize',fsize)
        title('Class Data','FontSize',fsize)
        
        figure
        if(length(unique(Class_Estimate)) == 2)
            gclrz = 'br';
        else
            gclrz = 'gbr';
        end
        %gscatter(Data_All(:,1),Data_All(:,2),Class_Estimate,gclrz);
        strlab = {'Class 1';'Class 2'};
        gscatter(Data_All(:,1),Data_All(:,2),strlab(Labels),'br');
        for cc = 1:length(cslist)
            hold on
            Surface3D(DataClass{cc}(:,1),DataClass{cc}(:,2),Difference{cc}, 'mesh');
        end
        hold off
        str = sprintf('Seperability (%f)',mean(Diff_All) );
        title(str,'FontSize',fsize)
        xlabel('State 1','FontSize',fsize)
        ylabel('State 2','FontSize',fsize)
        zlabel('W_{rbf}','FontSize',fsize,'interpreter','latex')
    elseif(SS == 3)
        figure
        for cc = 1:length(cslist)
            scatter3(DataClass{cc}(:,1),DataClass{cc}(:,2),DataClass{cc}(:,3),plotsize3,Difference{cc});
            hold on
        end
        str = sprintf('Seperability w/ Discriminant (%f)',mean(AccuracySeperable) );
        title(str)
        xlabel('x1')
        ylabel('x2')
        zlabel('x3')
        colormap cool
        colorbar
        figure
        if(length(unique(Class_Estimate)) == 2)
            gclrz = 'br';
        else
            gclrz = 'gbr';
        end
        gscatter3(Data_All(:,1),Data_All(:,2),Data_All(:,3),Class_Estimate,gclrz);
        title('class estimate')
        xlabel('x1')
        ylabel('x2')
        zlabel('x3')
    else
        warning('cannot plot this many dimensions');
        figure
        for cc = 1:length(cslist)
            scatter3(DataClass{cc}(:,1),DataClass{cc}(:,2),DataClass{cc}(:,3),plotsize3,Difference{cc});
            hold on
        end
        str = sprintf('Seperability w/ Discriminant (%f)',mean(AccuracySeperable) );
        title(str)
        xlabel('x1')
        ylabel('x2')
        zlabel('x3')
        colormap cool
        colorbar
        figure
        if(length(unique(Class_Estimate)) == 2)
            gclrz = 'br';
        else
            gclrz = 'gbr';
        end
        gscatter3(Data_All(:,1),Data_All(:,2),Data_All(:,3),Class_Estimate,gclrz);
        title('class estimate')
        xlabel('x1')
        ylabel('x2')
        zlabel('x3')
    end
    
end


%find good data where relative RBF is high enough
good_data = Data_All(Diff_All > thresh,:);
good_data_class = Class_All(Diff_All > thresh,:);
good_data_Diff = Diff_All(Diff_All > thresh,:);

%Fit linear or gaussian mixtures from seperable data
Model = [];
allthresh = [];

for cc = 1:length(cslist)
    Dfit = good_data(good_data_class == cslist(cc),:);
    DiffFit = good_data_Diff(good_data_class == cslist(cc),:);

    %T = clusterdata(Dfit,'maxclust',3);
    idc = 0;
    if(length(Dfit) > 10)
        eva = evalclusters(Dfit,'kmeans','gap','KList',[1:kmeanlimit]);
        km(cc) = eva.OptimalK; %assume this is solved
        %T = NormalizedGraphCuts(Dfit,groupings(cc));
        [meansk,T] = kMeansIterative(Dfit,km(cc));
        for ii = 1:max(T)
            Dclust = Dfit(T==ii,:);
            DiffClust = DiffFit(T==ii,:);
            if(length(Dclust) > 10)
                idc = idc + 1;
                %pd = fitdist(Dclust,'Normal');
                Model{cc}.cluster{idc}.mean = mean(Dclust);
                Model{cc}.cluster{idc}.sigma = cov(Dclust);
                Model{cc}.cluster{idc}.scale = gaussianScale(Model{cc}.cluster{idc}); %1/gaussianProbMV(Model{cc}.cluster{idc}.mean,Model{cc}.cluster{idc}.sigma,Model{cc}.cluster{idc}.mean);
                Model{cc}.cluster{idc}.clusterdata = Dclust;
                %compute probabilities for cluster data
                Model{cc}.cluster{idc}.prob = gaussianEval(Dclust,Model{cc}.cluster{idc});
                
                %find 95% for threshold
%                 [DiffMin,DI] = sort(abs(DiffClust - thresh));
%                 nd = round(length(DiffMin)*0.95);
%                 sortprob = Model{cc}.cluster{idc}.prob(DI,:);
%                 allthresh = [allthresh; mean(sortprob(1:nd))];
                
                guesser = 2*min(Model{cc}.cluster{idc}.prob); %0.1;
                allthresh = [allthresh; guesser];
            end
        end
    else
        warning('class %d does not contain any seperable data', cc)
    end
    Model{cc}.TotalClusters = idc;
    Model{cc}.SeperableData = Dfit;
    Model{cc}.Difference = good_data_Diff(good_data_class == cslist(cc),1);
    Model{cc}.ClassLabel = cslookup{cc};
    Model{cc}.ThresholdSep = mean(allthresh)/2;
end




%plot if need be
if(strcmp(ploton,'ploton'))
    if(SS == 1)
        figure
        gscatter(good_data,zeros(length(good_data),1),good_data_class,'br');
        for cc = 1:length(cslist)
            for idc = 1:Model{cc}.TotalClusters
                hold on
                scatter(Model{cc}.cluster{idc}.clusterdata,Model{cc}.cluster{idc}.prob,'g.')
            end
        end
        hold off
        title('Seperable Data and new Fit')
        xlabel('data')
        ylabel('Probability')
    elseif(SS == 2)
        figure
        gscatter(good_data(:,1),good_data(:,2),good_data_class,'br');
        for cc = 1:length(cslist)
            for idc = 1:Model{cc}.TotalClusters
                hold on
                Surface3D(Model{cc}.cluster{idc}.clusterdata(:,1),Model{cc}.cluster{idc}.clusterdata(:,2),Model{cc}.cluster{idc}.prob);
            end
        end
        hold off
        title('Seperable Data and new Fit')
        xlabel('x1')
        ylabel('x2')
        zlabel('probability')
    elseif(SS == 3)
        figure
        for cc = 1:length(cslist)
            for idc = 1:Model{cc}.TotalClusters
                scatter3(Model{cc}.cluster{idc}.clusterdata(:,1),Model{cc}.cluster{idc}.clusterdata(:,2),Model{cc}.cluster{idc}.clusterdata(:,3),plotsize3,Model{cc}.cluster{idc}.prob);
                hold on
            end
        end
        hold off
        title('Seperable Data and new new fit (color)')
        xlabel('x1')
        ylabel('x2')
        zlabel('x3')
        colormap cool
        colorbar
    else
        warning('cannot plot this many dimensions, only plotting 3')
        figure
        for cc = 1:length(cslist)
            for idc = 1:Model{cc}.TotalClusters
                scatter3(Model{cc}.cluster{idc}.clusterdata(:,1),Model{cc}.cluster{idc}.clusterdata(:,2),Model{cc}.cluster{idc}.clusterdata(:,3),plotsize3,Model{cc}.cluster{idc}.prob);
                hold on
            end
        end
        hold off
        title('Seperable Data and new new fit (color)')
        xlabel('x1')
        ylabel('x2')
        zlabel('x3')
        colormap cool
        colorbar
    end
    
end

end