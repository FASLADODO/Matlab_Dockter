function [PL,Y_Est] = LinearGaussianOnline(Data_Online,Y_Online,Model)
%http://www.cs.cmu.edu/~guestrin/Class/10708/slides/gaussians-kf.pdf

    %Now get mean from data with params
    Y_Est = Data_Online*Model.paramLG;
    
    %checks if in range (since least squares goes on forever)
    mask =  inRange(Data_Online,Model.min,Model.max);
    
    %get linear gasussian for online data
    PL = (1/ ( (2*pi)^(Model.ss/2) .*sqrt(norm(Model.sigL)) ) ).* exp(-(1/2).* QuadraticForm(Y_Online - Y_Est,Model.sigin));
%     PL = sqrt(norm(Model.sigL) / (2*pi) ) .* exp(-(1/2).* QuadraticForm(Y_Online - Y_Est,Model.sigin));
    
    PL = PL';
    %if not in range set to zero
    PL = double(mask).*PL;
    
end