function X = subSampleVec(X_in,sub_percent)

    %check format
    [n_o,dim_o]=size(X_in);
    if(dim_o > n_o) %check on format should be (3 columns x n rows)
        X_in = X_in';%resize
        [n_o,dim_o]=size(X_in);
    end

    %subsample original points
    ss = round(n_o*sub_percent);
    
    %get step size
    steps = 1/sub_percent;
    subsam = round(1:steps:n_o)
    
    X = X_in(subsam,:);
    
    [n,dim]=size(X);
    
    disp(sprintf('Using %i , out of %i data points \n', n, n_o ));


end