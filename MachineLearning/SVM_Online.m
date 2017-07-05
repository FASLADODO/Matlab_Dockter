function [Accuracy,Classify] = SVM_Online(X,Y,W,B)
    %Test SVM accuracy
    Classify = [];
    %test classify
    for i = 1:length(X)
       Z = X(i,:)';
       sumz = W*Z + B;
       class = sign(sumz);

       Classify = [Classify; Y(i), class];

    end

    correct = double(Classify(:,1) == Classify(:,2));
    Accuracy = sum(correct)/length(correct);

end