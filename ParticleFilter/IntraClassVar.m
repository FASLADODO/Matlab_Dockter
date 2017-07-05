function [ var ] = IntraClassVar( class1, class2 )
    combined = (class1 + class2)/2;
    var = 0;
    avg = 0;
    %Strange average computation from pdf
    for ii = 1:length(combined)
        avg = avg + ii*combined(ii);
    end
    %compute variance of vector
    for ii = 1:length(combined)
        var = var + combined(ii)*((ii-avg)^2);
    end
end


