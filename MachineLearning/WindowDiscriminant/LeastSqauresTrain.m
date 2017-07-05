function params = LeastSqauresTrain(Xon,Yon,labels)
    
    classes = unique(labels);

    for kk = 1:length(classes)
        %Get seperate data for each class
        datac{kk} = Xon(labels == classes(kk),:);
        yc{kk} = Yon(labels == classes(kk),:);
        
        params{kk} = pinv(datac{kk})*yc{kk};
    end
    
end