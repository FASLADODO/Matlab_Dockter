function [D,alpha] = AdaBoostUpdate(D,Y,H_x)
    %update weights for weak learners
    %Taken from page 122 Mohri Fundamentals of Machine Learning
    %D: previous iteration weights
    %error: classification error from previous learner
    %Y: true class label
    %H_x: Class estimate label
    
    %classification error
    error = ClassificationError(Y,H_x);
    
    %update weights
    if(error > eps)
        alpha = 0.5 * log( (1-error)/error ); % shift weight
    else
        alpha = 0.5 * log( 1/eps );
    end
    
    D = D.*exp( -alpha*Y.*H_x ) ;%/ Zn;
    
    %Normalization
    %Zn = 2*( error*(1-error) )^(1/2); 
    D = D./sum(D);

end