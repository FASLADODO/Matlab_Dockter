function [PL,meany] = LinearGaussian(Data_Train,Y_Train,Data_Online,Y_Online)

    %LS Parameters
    params = pinv(Data_Train)*Y_Train
    meany = Data_Train*params;
    
    %Covariance from mean
    sigL = cov(Y_Train-meany);
    sigin = inv(sigL); %inverse
    ss = length(sigL); %Covariance size
    
    %Now get mean from data with params
    meanyo = Data_Online*params;
    
    
    %get linear gasussian for online data
    PL = (1/ ( (2*pi)^(ss/2) .*sqrt(norm(sigL)) ) ).* exp(-(1/2).* QuadraticForm(Y_Online - meanyo,sigin));
    PL = PL';

end