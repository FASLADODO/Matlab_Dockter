function classify = DPR_Online(PDFS,X)

classes = length(PDFS);
for cc = 1:classes
    points(cc) = length(PDFS{cc}.lin);
end
order = length(X(1,:));

NN = length(X);

%loop through all points
for ii = 1:NN
    probs = [];
    %loop through all classes
    for cc = 1:classes
        idx = [];
        %loop through all orders
        for oo = 1:order
            %get index
            idx(oo) = round( (X(ii,oo) - PDFS{cc}.origin(oo) )/ PDFS{cc}.spacing(oo) ) + 1;
            if(idx(oo) < 1)
                idx(oo) = 1;
            end
            if(idx(oo) > points(cc))
               idx(oo) = points(cc); 
            end
        end
        %get the probability at the current location
        B_ind = num2cell(idx);
        probs(cc) = PDFS{cc}.prob_grid(B_ind{:});
    end
    %get class from ratios
    classify(ii) = pairwise_ratio(probs,0);
    
end


end