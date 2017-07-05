function [ var ] = Variance( probs )
%get variance from probability vector
    var = 0;
    avg = 0;
    %Strange average computation from pdf
    for ii = 1:length(probs)
        avg = avg + ii*probs(ii);
    end
    %compute variance of vector
    for ii = 1:length(probs)
        var = var + probs(ii)*((ii-avg)^2);
    end

end

