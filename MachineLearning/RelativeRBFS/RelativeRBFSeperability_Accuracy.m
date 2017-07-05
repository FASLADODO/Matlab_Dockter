function [Variations,BestSep,Max_Diff_All] = RelativeRBFSeperability_Accuracy(Data,Labels,testcolumns,columnlabels,minchoosek,ploton)
%Use this to check the max / mean seperability for a given state space
%Using classification accuracy as measure

%Data = class data to train on (rows are samples, columns are dimensions)
%Labels = column of unique class lables with same length as Data [1;2...]

%ploton = 'ploton' or 'plotoff'

[NN,SS] = size(Data);

%get all unique classes
if(isnumeric(Labels))
    cslist = unique(Labels);
    cslookup = strread(num2str(cslist'),'%s');
else
    cslookup = unique(Labels);
    cslist = 1:length(cslookup);
end

%get the multi combo nchoosek
maxk = length(testcolumns);
kmulti = minchoosek:maxk;
nkcombos = nchoosekmulti(testcolumns,kmulti);

plotsize3 = 10;

%for finding the bestest
oldmax = 0;
BestSep.bestcolumns = [];
BestSep.bestnumcol = 1;
BestSep.bestcollabels = {};
BestSep.max = 0;
BestSep.Difference = [];
BestSep.DataClass = [];

%All Data and Diff
Data_All = [];
Diff_All = [];
Class_All = [];

Max_Diff_All = [];

mr = max(cslist) + 1;

Variations = [];
%loop through all possible column/state combinations
for nk = 1:length(nkcombos)
    
    %loop through all combos of given size
    for ii = 1:size(nkcombos{nk},1)
        columni = nkcombos{nk}(ii,:)
        
        %loop through all classes
        Data_All_Temp = [];
        Diff_All_Temp = [];
        Class_All_Temp = [];
        Max_Diff = [];
        Mean_Diff = [];
        
        %Get a test threshold for this combo
        SNK = length(columni);
        %thresh = log10(2)/(SNK); %SNK^1.6
        thresh = 0.001; %2*10^(-((pi/2)*SNK));

        %Find all relative RBFS for class data
        for cc = 1:length(cslist)
            %current class and opposite class
            DONt = Data(Labels == cslist(cc),:);
            DOFFt = Data(Labels ~= cslist(cc),:);
            DON = DONt(:,columni);
            DOFF = DOFFt(:,columni);
            
            %optimal bandwidth (rule of thumb)
            sig = norm(std(DON));
            nn = size(DON,1);
            bw = 1.06*sig*(nn^(-1/5));

            %compute dem rbfs
            [pwithin,pbetween] = ComputeDiscriminateRBF(DON,DOFF,bw);

            Difference{cc} = ComputeRBFDifference(pwithin,pbetween);
            %fo lata
            SavePR{cc}{cc} = Difference{cc};
            SavePR{cc}{mr-cc} = ComputeRBFDifference(pbetween,pwithin);
            DataClass{cc} = DON;

            %extra stuff
            Data_All_Temp = [Data_All_Temp; DON];
            Diff_All_Temp = [ Diff_All_Temp; Difference{cc}];
            Class_All_Temp = [Class_All_Temp; ones(nn,1)*cc ];
            
            %Store everything
            Variations.nksize{nk}.row{ii}.class{cc}.data = DON;
            Variations.nksize{nk}.row{ii}.class{cc}.Diff = Difference{cc};
            
        end
        
        %TRY THIS!
        thresh = 0.3*max(Diff_All_Temp);
        
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

        metric =  mean(AccuracySeperable)
        if(metric > oldmax)
            oldmax = metric;
            BestSep.bestcolumns = columni;
            BestSep.bestnumcol = nk;
            BestSep.bestcollabels = columnlabels(columni);
            BestSep.max = Max_Diff;
            BestSep.Difference = Difference;
            BestSep.DataClass = DataClass;
            Data_All = Data_All_Temp;
            Diff_All = Diff_All_Temp;
            Class_All = Class_All_Temp;
        end
    end
end


%function for pretty plots
RBFNicePlots(BestSep.DataClass,BestSep.Difference,Data_All,Class_All,BestSep.bestcollabels);




end