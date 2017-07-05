function [P,Lm] = LogRegOnline(Data,Model)
%LogRegOnline: compute class probabilities online from Logisitic Regression
%Model comes from LogRegTrain()

%P: probability matrix (each column corresponds to P for each class
%Lm: is overall logisitic mapping

    mapping = Data*Model.Params; %map to logistic space
    Lm = Logistic(mapping);

    if(size(Lm,2) > 1)
        P = Lm;
    else
        %get probabilities from each class
        P1 = 1-Lm;
        P2 = Lm;
        %return them in columns
        P = [P1,P2];
    end
end