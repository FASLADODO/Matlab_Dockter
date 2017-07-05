function [x_hat_m, P_m] = SLAM_EKF_Propogate(x_hat_p, P_p, v_m, sigma_v, w_m, sigma_w, dt, Nc)
%Rod Dockter
%Code for EKF propogation
%CSCI 5552
%Homework 4 problem 2
%Ported from c++ code for project


% x_hat(k+1|k) = f[x_hat(k|k),u(k),0]
% P(k+1|k) = Phi(k)*P(k|k)*Phi(k)' + G(k)*Q(k)*G(k)'
x_hat_m = x_hat_p;

% Define the linearized state transition matrix
Phi = [1 0 -v_m*dt*sin(x_hat_p(3,1));...
    0 1 v_m*dt*cos(x_hat_p(3,1));...
    0 0 1];

% Define the linearized noise matrix
G = [-dt*cos(x_hat_p(3,1)) 0;...
    -dt*sin(x_hat_p(3,1)) 0;...
    0 -dt];

if Nc ~= 0
    Phi = [Phi, zeros(3,2*Nc);...
        zeros(2*Nc,3), eye(2*Nc)];
    G = [G;...
        zeros(2*Nc,2)];
end

% Define the system noise coveriance matrix
Q = (v_m)^2*[sigma_v^2 0;...
            0 sigma_w^2];

% Propogate the trajectory using the previous state estimate and the
% measured velocity
x_hat_m(1,1) = x_hat_p(1,1) + v_m*dt*cos(x_hat_p(3,1));
x_hat_m(2,1) = x_hat_p(2,1) + v_m*dt*sin(x_hat_p(3,1));
x_hat_m(3,1) = x_hat_p(3,1) + w_m*dt;

% Propogate the covariance using the previous covariance estimate and noise
% covariance
P_m = Phi*P_p*Phi' + G*Q*G';

end
