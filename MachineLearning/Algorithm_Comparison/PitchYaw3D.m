function [alpha,beta] = PitchYaw3D(X)
%returns alpha and beta angle columns for 3D vectors
%X row wise data matrix, columns are dimensions

    [~,SS] = size(X);
    if(SS ~= 3)
       error('wrong number of columns in X') 
    end
    alpha = atan2(X(:,2),X(:,1));
    beta = atan2(X(:,3),X(:,1));

end