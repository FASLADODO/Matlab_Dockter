function E = ErrorSOS(Y,Ybar)
%Sum of squares error between true and estimate

    if(length(Y) ~= length(Ybar))
       error('vectors must be same length') 
    end
    
    E = (1/2)*sum( abs(Y-Ybar).^2 );
end