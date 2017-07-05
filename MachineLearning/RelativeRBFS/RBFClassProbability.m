function [prob] = RBFClassProbability(TrainData,OnlineData)
    % Non-Class relative
    %simply compute the probability for each point in online data using
    %RBFs from training data
    
    [NT,ST] = size(TrainData);
    [NO,SO] = size(OnlineData);
    
    %optimal bandwidth (rule of thumb)
    sig = norm(std(TrainData));
    n = length(TrainData);
    bw = 1.06*sig*(n^(-1/5))

    %compute dem rbfs
    onlinep = rbfpdist2(OnlineData,TrainData,bw);
    prob = double(onlinep/NO);
end