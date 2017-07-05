function Pprime = Transform2x2(P,T)
    %Adds a ones column and maps data through a transformation matrix
    
    [NN,SS] = size(P);
    
    if SS ~= 2
       error('coordinate matrix must be n x 2') 
    end
    
    P = [P, ones(NN,1)];
    
    Pnew = T * P' ;
    
    Pprime = Pnew';

end