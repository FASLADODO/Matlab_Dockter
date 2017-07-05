function [Model,GoodData,thresh] = SimpleRelativeRBFProbs(Difference,ClassData,ProbData)

    %loop through all classes
    GoodData = [];
    Model = [];

    %Find all relative RBFS for class data
    for cc = 1:length(Difference)
        %current class
        DON = ClassData{cc};
        if(length(DON) > 0)
            Diff = Difference{cc};
            %thresh = 0.1*max(Diff);
            thresh(cc) = 0.2*prctile(Diff,95)

            %find all data above a threshold
            idx = find(Diff > thresh(cc));

            %stash it
            GoodData{cc} = ClassData{cc}(idx,:);

            %Make a model
            Model{cc}.mean = mean(GoodData{cc});
            Model{cc}.sigma = cov(GoodData{cc});
            Model{cc}.prob = mvnpdf(GoodData{cc},Model{cc}.mean,Model{cc}.sigma); 
            Model{cc}.scale = 1/ max(Model{cc}.prob);
            Model{cc}.thresh = thresh(cc);
        end
    end
    
    
end