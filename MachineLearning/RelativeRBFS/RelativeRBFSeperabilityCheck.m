function [Variations,BestSep,BestSepClass,Data_All,Diff_All] = RelativeRBFSeperabilityCheck(Data,Labels,testcolumns,columnlabels,minchoosek,ploton)
%Use this to check the max / mean seperability for a given state space

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

%All Data and Diff
Data_All = [];
Diff_All = [];
Class_All = [];

Max_Diff_All = [];

Variations = [];
%loop through all possible column/state combinations
for nk = 1:length(nkcombos)
    oldmax(nk) = 0; %do it for each number of states
    
    %loop through all combos of given size
    for ii = 1:size(nkcombos{nk},1)
        columni = nkcombos{nk}(ii,:);
        
        %loop through all classes
        Data_All_Temp = [];
        Diff_All_Temp = [];
        Class_All_Temp = [];
        Max_Diff = [];
        
        %Get a test threshold for this combo
        SNK = length(columni);
        %thresh = log10(2)/(SNK); %SNK^1.6
        thresh = 2; %2*10^(-((pi/2)*SNK));

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
            DataClass{cc} = DON;

            %extra stuff
            Data_All_Temp = [Data_All_Temp; DON];
            Diff_All_Temp = [ Diff_All_Temp; Difference{cc}];
            Class_All_Temp = [Class_All_Temp; ones(nn,1)*cc ];
            
            %Store everything
            Variations.nksize{nk}.row{ii}.class{cc}.data = DON;
            Variations.nksize{nk}.row{ii}.class{cc}.Diff = Difference{cc};
            
            %try only taking diffs above threshold
            %Max_Diff(cc) = mean(Difference{cc} > thresh);
            %Max_Diff(cc) = sum(Difference{cc} > thresh);
            %Max_Diff(cc) = mean(Difference{cc});
            %Max_Diff(cc) = mean(gooddiff);
            Max_Diff(cc) = mean(Difference{cc});
        end
        Max_Diff_All =[Max_Diff_All; Max_Diff,SNK];
        
        %Store everything
        Variations.nksize{nk}.row{ii}.columns = columni;
        Variations.nksize{nk}.row{ii}.ncol = nk;
        Variations.nksize{nk}.row{ii}.labels = columnlabels(columni);
        
        metric =  weightingScore(Max_Diff);
        if(metric > oldmax(nk))
            oldmax(nk) = metric;
            BestSepClass{nk}.bestcolumns = columni;
            BestSepClass{nk}.bestnumcol = nk;
            BestSepClass{nk}.bestcollabels = columnlabels(columni);
            BestSepClass{nk}.max = Max_Diff;
            BestSepClass{nk}.metric = metric;
            BestSepClass{nk}.Difference = Difference;
            BestSepClass{nk}.DataClass = DataClass;
            Data_All{nk} = Data_All_Temp;
            Diff_All{nk} = Diff_All_Temp;
            Class_All{nk} = Class_All_Temp;
        end
    end
end

%also give back the best of all states
[~,bestsepid] = max(oldmax);
BestSep = BestSepClass{bestsepid};

figure
scatter3(Max_Diff_All(:,1),Max_Diff_All(:,2),Max_Diff_All(:,end))
%axis([0,1,0,1,0,5])
xlabel('Seperability 1')
ylabel('Seperability 2')
zlabel('# States')

%function for pretty plots
%RBFNicePlots(BestSep.DataClass,BestSep.Difference,Data_All,Class_All,BestSep.bestcollabels);




end