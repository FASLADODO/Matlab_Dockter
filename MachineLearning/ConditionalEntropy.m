function H_Y_X = ConditionalEntropy(Y,X,threshold)
% Rods implementation
%compute conditional entropy of Y given X
%"Entropy of Y given X"
%X can be real valued or categorical
%if X is categorical threshold should not be specified
%if X is numerical, then threshold will split X into 0/1
%Y = class labels
%optional threshold to classify X
%slide 9-10
%http://www.cs.cmu.edu/~guestrin/Class/10701-S06/Handouts/recitations/recitation-decision_trees-adaboost-02-09-2006.ppt

    %check sizes
    [~,SX] = size(X);
    if(SX > 1)
       warning('ignoring additional columns of X') ;
       X = X(:,1);
    end
    
    %in case there are more than two class labels
    classes = unique(Y);
    
    %choose mode
    %type is normally 'categorical';
    if(nargin == 3)
       %if 'numerical' use threshold to classify into 0,1
       temp = X > threshold;
       X = temp; %Now X is 0 or 1
    end

    %now lets compute conditional entropy
    cats = unique(X); %possible categories ([0,1] for numerical)
    H_X = [];
    P_X = [];
    for ss = 1:length(cats) %loop through all values in X (categories)
       X_C = X == cats(ss); %0 or 1 for given category
       P_X(ss) = mean(X_C); % probability of given category
       Y_X = Y(X_C); %class labels in given category
       
       %compute probability of Y among only those records which match X
       p = [];
       for cc = 1:length(classes) %loop through all classes
           temp = Y_X == classes(cc); % indices for each class given category
           p(cc) = mean(temp); %probability of given class for given category
       end
       H_X(ss) = Entropy(p); %Entropy of class membership given X
    end
    H_Y_X = sum(P_X.*H_X); %compute conditional entropy
end