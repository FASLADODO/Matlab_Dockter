function [within,between] = ComputeDiscriminateRBF(Dwithin,Dbetween,bandwidth)
%compute inter and intra RBFs for discriminate analysis
%See ProbabilityRadialBasis.m
%Dwithin is all data within current class (rows are samples)
%Dbetween is data from other class (rows are samples)
%bandwidth is the RBF scaling bandwidth

    %scale by lengths
    nwithin = size(Dwithin,1);
    nbetween = size(Dbetween,1);

    %compute optimal bandwidth if not supplied
    if(nargin == 2)
        Data = [Dwithin; Dbetween];
        nb = sum([nwithin,nbetween]);
        sig = norm(std(Data));
        bandwidth = 1.06*sig*(nb^(-1/5)) ;
    end

    %get pdists
    %get within and between class
    within_us = rbfpdist2(Dwithin,Dwithin,bandwidth);
    between_us = rbfpdist2(Dwithin,Dbetween,bandwidth);
    
    %get scale from integral estimate
    %scalew = ScaleRBF(Dwithin,within_us);
    %scaleb = ScaleRBF(Dbetween,between_within_us);
    %[scalew,~,~] = integralRBF(Dwithin,bandwidth);
    %[scaleb,~,~] = integralRBF(Dbetween,bandwidth);
    scalew = nwithin;
    scaleb = nbetween;
%     scalew = max(within_us);
%     scaleb = max(between_us);
    
    %now scale the within and between estimates
    within = double(within_us/scalew);
    between = double(between_us/scaleb);
end