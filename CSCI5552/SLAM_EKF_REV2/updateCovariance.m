function [P] = updateCovariance(P,K,H,S)
%tame the covariance

%innovation covariance
Z = H*P*H' + S;

%adjust covariance
P = P - K*Z*K';

end