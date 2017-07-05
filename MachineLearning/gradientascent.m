function theta = gradientascent(theta, alpha, Y, X)
    %simple woodruff hoff learning
    %alpha is learning rate 0-1
    
    m = length(Y);
    thetaLen = length(theta);
    
    %X*theta = x_0*th_0 + x_1*th_1 + ... + x_n*th_n
    for k = 1:m
        h_th = sigmoid(X(k,:),theta);
        temp = (X*theta - h_th);

        for i=1:thetaLen
            tempVal(i,1) = sum(temp.*X(k,i));
        end

        theta = theta + alpha*tempVal;
    end

end