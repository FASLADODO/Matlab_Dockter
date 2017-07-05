function [K, gotime] = computeKalmanGain(P,H,S)
    %compute kalman gain given jacobians, covariance and landmark updates
    %tells us how much we trust odom and landmarks

    %innovation covariance
    Z = H*P*H' + S;
    
    gotime = 0;
    if(cond(S) < 80) %to make sure the inverse is well behaved
        %Kalman Gain
        K = P*H'*inv(Z);
        gotime = 1;
    else
        K = 0;
    end
    
end