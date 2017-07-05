function E = ErrorRMS(Y,Ybar)
    
    if(length(Y) ~= length(Ybar))
       error('vectors must be same length') 
    end
    [NN,~] = size(Y);
    
    E = sqrt( (1/NN)*sum( abs(Y-Ybar).^2 ) );

end