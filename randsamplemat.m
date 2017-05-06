function sample = randsamplemat(mat,percent)
%randsamplemat(mat,percent): helper function to get a percentage of random
%samples from an existing mxn matrix, where m rows is the longer dimension.
%Returns the randomly sampled matrix

    %get length
    vec_length = length(mat);
    
    %get number of points to sample
    kk = round(percent*vec_length);
    
    %get sub sample indices
    y = randsample(vec_length,kk);
    
    %return sub sampled vector
    sample = mat(y,:);

end