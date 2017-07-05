function XP = MapValues(X,oldvalues,newvalues)
%MapValues: takes all values given in 'oldvalues' and maps them to those
%in 'newvalues' within matrix X

    if(length(oldvalues) ~= length(newvalues) )
       error('mapping vectors need to be the same length'); 
    end

    XP = zeros(size(X));
    for ii = 1:length(oldvalues)
        if(iscell(X))
            XP( strcmp(X,oldvalues(ii)) ) = newvalues(ii);
        else
            XP(X == oldvalues(ii)) = newvalues(ii);
        end
    end

end