function [PC] = GaussianClassifierOnline(Data,Model)
%Taken from CSCI 5525 Leg04-gendisc

    [NN,SS] = size(Data);
    PC = sigmoidFunction(Data * Model.W + repmat(Model.w0,NN,1) );
end