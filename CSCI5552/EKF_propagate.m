function [x_hat_min, P_min] = EKF_propagate(x_hat_plus, P_plus, v_m, w_m, Q, dt)
%%Rod Dockter
%%Homework 3
%%problem 3
%%csci 5552
%%This is just the standard EKF propogate equations using the sample code
%%from homework 2
%%This gets called by thetwoD_EKF main program

% state propagation
x_hat_min = x_hat_plus + dt*[v_m*cos(x_hat_plus(3));  
                             v_m*sin(x_hat_plus(3));
                             w_m];

% Jacobians
Phi = [ 1   0   -dt*v_m*sin(x_hat_plus(3));
        0   1    dt*v_m*cos(x_hat_plus(3));
        0   0    1                          ];


G = -dt*[cos(x_hat_plus(3)) 0;
         sin(x_hat_plus(3)) 0;
         0                  1];


% covariance propagation         
P_min = Phi*P_plus*Phi' + G*Q*G';