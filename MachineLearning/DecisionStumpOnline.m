function Class = DecisionStumpOnline(X,Thresh,Direction)
%classify in one dimension given a threshold

    [NN,SS] = size(X);

    if(SS > 1)
       warning('X has too many columns, using only the first')
       X = X(:,1);
    end

    %check thresh
    isabove = X >= Thresh;
    
    %see which class that is
    Class = Direction(isabove+1);
    
    %make it gooder
    if(~iscolumn(Class))
       Class = Class'; 
    end
end