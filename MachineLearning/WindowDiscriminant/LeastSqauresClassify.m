function Accuracy = LeastSqauresClassify(Xon,Yon,labels,params)
    
    classes = unique(labels);

    for kk = 1:length(classes)
        %Get seperate data for each class
        datac{kk} = Xon(labels == classes(kk),:);
        yc{kk} = Yon(labels == classes(kk),:);
        
        if(nargin == 3)
            params{kk} = pinv(datac{kk})*yc{kk};
        end
    end
    
    Ratios = [];
    Classified = [];
    for ii = 1:length(classes)
        for jj = 1:length(classes)
            errortemp = abs(yc{ii} - datac{ii}*params{jj});
            Errors{ii}(:,jj) = errortemp;
        end
        
        
        [~,idxclass] = min(Errors{ii},[],2);
        Classified = [Classified; idxclass];
    end
    
    Correct = Classified == labels;
    
    %Accuracy = sum(Correct)/length(Correct);

    ratios1 = log(Errors{1}(:,1) ./ Errors{1}(:,2));
    ratios2 = log(Errors{2}(:,1) ./ Errors{2}(:,2));
    correct1 = ratios1 <= 0;
    correct2 = ratios2 >= 0;
    c1acc = sum(correct1)/length(correct1)
    c2acc = sum(correct2)/length(correct2)
    
    Accuracy = [-ratios1, -ratios2];
end