function [cm] = CumMean(Data)
    %Data is a 1D column
    %cm is a running cumulative mean
    
    [NN,~] = size(Data);
    
    cm = cumsum(Data)./([1:NN]');
end