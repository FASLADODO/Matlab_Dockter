function KL = RelativeEntropy(p1,p2)

    if(length(p1) ~= length(p2))
       error('vectors must be same length'); 
    end
    KL = sum( p1 .* log( p1./(p2+eps)));
end