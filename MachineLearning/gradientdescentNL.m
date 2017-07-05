function params = gradientdescentNL(nl_func, params, Y, X, alpha, epochs, option, parametermapping)
    %simple woodruff hoff learning
    %alpha is learning rate 0-1
    %epochs is number of training epochs
    %X is matrix, columns are states, rows are samples
    %X must have the same number of columns as length(params) (repeat if
    %neccesary
    %Y is a column vector
    %params are 1D column array [theta1;theta2;theta3,...]
    %define function as:
    %nl_func = @(X,Param) Param(1)*sin(X(:,1)) + X(:,2) * Param(2);
    %option = 'batch' or 'stochastic'
    %parametermapping which columns affect which parameters (optional)
    
    m = length(Y);
    thetaLen = length(params);
    
    %default args
    if(nargin == 6)
        parametermapping = 1:thetalen;
    end
    if(nargin == 6)
        option = 'batch';
    end
    
    for tt = 1:epochs
        if(option == 'batch')
            %update theta once for all samples

            %compute estimate of output
            %X*theta = x_0*th_0 + x_1*th_1 + ... + x_n*th_n
            Y_bar = nl_func(X,params);

            %get residual error
            temp = (Y_bar - Y);

            %get error times each state
            for i=1:thetaLen
                errorScale(i,1) = sum(temp.*X(:,parametermapping(i)));
            end

            %update parameters
            params = params - (alpha/m)*errorScale;

        elseif(option == 'stochastic')
            %loop through and update theta for each sample
            for k = 1:m

                Y_bar = nl_func(X(k,:),params);
                temp = Y_bar - Y(k,1);


                for i=1:thetaLen
                   errorScale(i,1) = sum(temp.*X(k,parametermapping(i)));
                end

                params = params - alpha*errorScale;
            end
        end
    end
end