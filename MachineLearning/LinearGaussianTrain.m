function [Model] = LinearGaussianTrain(Data_Train,Y_Train)

% Data_Train: samples are rows, dimensions are columns
% Y_train: output variable (same length as Data_Train)

%See testLinearGaussian for example
%http://www.cs.cmu.edu/~guestrin/Class/10708/slides/gaussians-kf.pdf

    %LS Parameters
    Model.paramLG = pinv(Data_Train)*Y_Train;
    Model.meanLS = Data_Train*Model.paramLG;
    
    %Covariance from mean
    Model.sigL = cov(Y_Train-Model.meanLS);
    Model.sigin = inv(Model.sigL); %inverse
    Model.ss = length(Model.sigL); %Covariance size
    
    %Limits (may come in handy)
    Model.min = min(Data_Train,[],1);
    Model.max = max(Data_Train,[],1);
end