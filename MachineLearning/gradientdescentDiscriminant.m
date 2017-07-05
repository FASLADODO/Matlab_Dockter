function [theta1, theta2] = gradientdescentDiscriminant(theta1, theta2, alpha, Y1, X1, Y2, X2)
    %learn parameters using gradient descent for correct class, gradient
    %ascent for incorrect
    
    %SOMETHING LIKE THIS< MAYBE NEEDS TO CHANGE
    
    m = length(Y1);
    thetaLen = length(theta1);
    
  %update theta once for all samples
    %X*theta = x_0*th_0 + x_1*th_1 + ... + x_n*th_n
    temp1_1 = (X1*theta1 - Y1);
    temp1_2 = (X2*theta1 - Y2);
    temp2_1 = (X1*theta2 - Y1);
    temp2_2 = (X2*theta2 - Y2);

    for i=1:thetaLen
        tempVal1(i,1) = sum(temp1_1.*X1(:,i)) - sum(temp1_2.*X1(:,i));
        tempVal2(i,1) = sum(temp2_2.*X2(:,i)) - sum(temp2_1.*X2(:,i));
    end

    theta1 = theta1 - (alpha/m)*tempVal1;
    theta2 = theta2 - (alpha/m)*tempVal2;
end