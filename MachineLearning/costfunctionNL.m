function J = costfunctionNL(Y,X,nl_func,params)
    Y_bar = nl_func(X,params);
    J = (1/2) * sum( (Y - Y_bar).^2 );
end