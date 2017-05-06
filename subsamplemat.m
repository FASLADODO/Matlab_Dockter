function sample = subsamplemat(mat,percent)
%subsamplemat(mat,percent): helper function to get a sub sample
%from an existing mxn matrix, where m rows is the longer dimension
%Returns the sampled matrix

    %get length
    vec_length = length(mat);
    
    %get number of points to sample
    kk = round( vec_length / (percent*vec_length) );
    
    %return sub sampled vector
    sample = mat(1:kk:end,:);

end