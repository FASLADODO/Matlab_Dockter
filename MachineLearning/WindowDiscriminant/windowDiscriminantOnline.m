function LLR = windowDiscriminantOnline(DataOn, Models, classes, classthresh)
    %DataOn = X1 (sample set not training);
    %Models = struct of gaussian mixture model parameters for each class
    %classes = [1,2,3..] array of unique classes
    %classthresh = classification threshold below which probability is 0.5
    
    [NN,SS] = size(DataOn);
   
    Pall = [];

    for ii = 1:length(DataOn)
        datatemp = DataOn(ii,:);
        Pt = [];
        for cc = 1:length(classes)
            for mm = 1:length(Models{cc}.cluster)
                %%get scaled probability for each model at each data point
                Pt(mm,cc) = Models{cc}.cluster{mm}.scale * gaussianProbMV(datatemp,Models{cc}.cluster{mm}.sigma,Models{cc}.cluster{mm}.mu);
            end
        end
        valt(ii,1) = max(Pt(:));
        [rowm,colm] = find(Pt == valt(ii));
        if(valt(ii)) >= classthresh
            Pall(ii,colm) = valt(ii);
            Pall(ii,AllBut(classes,colm)) = 1 - valt(ii);
        else
            Pall(ii,:) = [0.5,0.5];
        end
        
        if(ii > 2*SS) %once we have enough data
            [LLR(ii),combos] = LLRmulti(Pall);
        else
            LLR(ii) = 0; %we just dont know yet
        end
    end


end