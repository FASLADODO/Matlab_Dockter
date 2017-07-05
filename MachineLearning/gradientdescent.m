function theta = gradientdescent(theta, alpha, Y, X, option)
    %simple woodruff hoff learning
    %alpha is learning rate 0-1
    
    m = length(Y);
    thetaLen = length(theta);
    
    if(option == 'batch')
        %update theta once for all samples
        %X*theta = x_0*th_0 + x_1*th_1 + ... + x_n*th_n
        temp = (X*theta - Y);
        
        for i=1:thetaLen
            tempVal(i,1) = sum(temp.*X(:,i));
        end

        theta = theta - (alpha/m)*tempVal;
    
    elseif(option == 'stochastic')
        %loop through and update theta for each sample
        for k = 1:m
            
            temp = X(k,:)*theta_in - Y(k,1);
            
            for i=1:thetaLen
               tempval(i,1) = sum(temp.*X(k,i));
            end
            
            theta = theta - alpha*tempVal;
        end
    end
end