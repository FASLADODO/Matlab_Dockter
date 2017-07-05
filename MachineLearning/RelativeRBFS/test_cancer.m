load ovariancancer

[Labels,key] = Category2Numeric(grp);

[RANK,WEIGHT] = relieff(obs,grp,20);

%rods fancy discriminant
%[DPWEIGHT,RANKDP] = DiscriminantPearson(obs,Labels);


KEEP = 10;
X = double(obs(:,RANK(1:KEEP)));


figure
gscatter3(X(:,1),X(:,2),X(:,3),Labels);

%% just test simple plots

Data = X(:,[1,2]);

[Difference,ClassData,ProbData] = SimpleRelativeRBFTrain(Data,Labels);


if(size(Data,2) == 2)
    figure
    gscatter(X(:,1),X(:,2),Labels);
    hold on
    cc = 1;
    %Surface3D(ClassData{cc}(:,1),ClassData{cc}(:,2),ProbData{cc}.prob(:,1));
    hold on
    cc = 2;
    Surface3D(ClassData{cc}(:,1),ClassData{cc}(:,2),ProbData{cc}.prob(:,1));
    hold off
    title('probs')


    figure
    gscatter(X(:,1),X(:,2),Labels);
    hold on
    cc = 1;
    Surface3D(ClassData{cc}(:,1),ClassData{cc}(:,2),Difference{cc});
    hold on
    cc = 2;
    Surface3D(ClassData{cc}(:,1),ClassData{cc}(:,2),Difference{cc});
    hold off
    title('seperability')
end

if(size(Data,2) == 3)
    figure
    gscatter3(X(:,1),X(:,2),X(:,3),Labels);
    title('classes')
    
    figure
    cc = 1;
    scatter3(ClassData{cc}(:,1),ClassData{cc}(:,2),ClassData{cc}(:,3),10,Difference{cc});
    hold on
    cc = 2;
    scatter3(ClassData{cc}(:,1),ClassData{cc}(:,2),ClassData{cc}(:,3),10,Difference{cc});
    hold off
    title('seperability')
    colormap cool
    colorbar
end



%%

collabz = {'1';'2';'3';'4';'5';'6';'7';'8';'9'};
testcolumns = 1:9;

Data = X;

[Variations,BestSep,BestSepClass,Data_All,Diff_All] = RelativeRBFSeperabilityCheck(Data,Labels,testcolumns,collabz,1,'ploton');
BestSep.bestcollabels
columnrbf = BestSep.bestcolumns;

for(ii = 1:length(BestSepClass) )
  str = sprintf('states: %d, metric: %f', ii, BestSepClass{ii}.metric );
  disp(str)
end

testX = Data_All{3};
testDiff = Diff_All{3};

if(size(testX,2) == 3)

    figure
    scatter3(testX(:,1),testX(:,2),testX(:,3),10,testDiff);
    title('seperability')
    colormap cool
    colorbar
end
