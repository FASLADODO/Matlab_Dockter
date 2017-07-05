function [PL,Uest] = NonLinearGaussian(Data,U,paramstruct,nlfunction)

    %LS Parameters
    params = paramstruct.params;
    
    %Covariance from mean
    sig = paramstruct.sigma;
    sigin = inv(sig); %inverse
    ss = length(sig); %Covariance size
    
    %Now get mean from data with params
    Uest = nlfunction(Data,params);
    
    
    %get linear gasussian for online data
    PL = (1/ ( (2*pi)^(ss/2) .*sqrt(norm(sig)) ) ).* exp(-(1/2).* QuadraticForm(U - Uest,sigin));
    PL = PL';

end