function Y_bar = GeometricMargin(X,Y,W,B)
    %GeometricMargin: compute the minimum margin to decision boundary
    % Y={+-1}, X=[x1i,x2i,x3i; ] W = theta, B = theta_0
    gamma = [];
    for ii = 1:length(Y)
        %Same as the stupid W^T*X form
        gamma(ii) = Y(ii)*( X(ii,:)*(W./norm(W)) + B/norm(W) );
    end
    
    Y_bar = min(gamma(ii));
    
end

