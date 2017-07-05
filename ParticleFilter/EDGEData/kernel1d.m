function k = kernel1d(u)
%simple kernel function
    if( abs(u) <= 1 )
        k = 0.5;
    else
        k = 0;
    end
end