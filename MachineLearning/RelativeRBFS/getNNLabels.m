function targets = getNNLabels(Labels)

    NN = length(Labels);
    cslist = unique(Labels);
    SS = length(cslist);
    
    targets = zeros(NN,SS);
    
    for ii = 1:length(cslist)
        targets(Labels == cslist(ii),ii) = 1;
    end

end