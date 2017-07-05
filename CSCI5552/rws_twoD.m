function [x_true, y_true, theta_true, v_m, w_m ] = rws_twoD(N, dt, x_true_1, y_true_1, theta_true_1, v_true, w_true, sigma_v, sigma_w)
% Problem 3 assignment 2
%%Rod Dockter
% real world simulation for 2D-localization

% pre-allocate space
x_true = zeros(1,N);
y_true = zeros(1,N);
theta_true = zeros(1,N);
v_m = zeros(1,N);
w_m = zeros(1,N);

% initialize
x_true(1,1) = x_true_1;
y_true(1,1) = y_true_1;
theta_true(1,1) = theta_true_1;

% true trajectory according to equations 5.7 - 5.9
for i = 2:N
    %%True equations as discussed by hand in 3
    theta_true(1,i) = theta_true(1,i-1) + w_true(1,i-1)*dt;
    x_true(1,i) = x_true(1,i-1) + v_true(1,i-1)*dt*cos(theta_true(1,i-1));
    y_true(1,i) = y_true(1,i-1) + v_true(1,i-1)*dt*sin(theta_true(1,i-1));
    
end

% Odometry measurements
%%including the velocity term makes the noise velocity dependent as
%%specified
v_m = v_true + sigma_v.*randn(1,N-1);
w_m = w_true + sigma_w.*randn(1,N-1);
%%(No gps measurement)
