function [PL,meany] = LinearGaussianSimple(Data,Y,Params)

    %LS Parameters
    meany = Data*Params;
    
    %Covariance from mean
    sigL = cov(Y-meany);
    sigin = inv(sigL); %inverse
    ss = length(sigL); %Covariance size
    
    %get linear gasussian for online data
    PL = (1/ ( (2*pi)^(ss/2) .*sqrt(norm(sigL)) ) ).* exp(-(1/2).* QuadraticForm(Y - meany,sigin));
    PL = PL';

end