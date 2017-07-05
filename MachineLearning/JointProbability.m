function PXY = JointProbability(Px, Py)
    [pxg,pyg] = meshgrid(Px,Py);
    
    PXY = pxg.*pyg;

end