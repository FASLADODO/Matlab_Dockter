function [H_X,H_Y,H_XY] = JointEntropy(Px,Py)
    %get joint probability table
    PXY = JointProbability(Px, Py);
    
    %get individual entropies
    H_X = Entropy(Px);
    H_Y = Entropy(Py);
    
    %conditional entropy
    H_XY = -sum(sum(PXY.*log2(PXY+eps)));
    
    
    %mutual Information
    I = -sum(sum(PXY.*log2( (Px.*Py) / (PXY+eps) )));

end