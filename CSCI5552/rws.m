function [x_true, v_m, w_m, distance] = rws(N, dt, x_true_1, v_true, w_true, Xldmk, Yldmk)
%%Rod Dockter
%%Homework 3
%%problem 3
%%csci 5552

%%real world simulation for distance only
%%Gets called in twoD_EKF in oder to get robot locations and velocity,
%%angular velocity and distance measurements, distorted by noise.

% preallocate space
x_true = zeros(3,N);

% initialize
x_true(:,1) = x_true_1;

%%%distance to landmarks
%%distance = sqrt(xlandmark-xposition^2+ylandmark-yposition^2)
dstlmk = zeros(3,N);
dstlmk(1,1)=sqrt(((Xldmk(1)-x_true_1(1,1))^2)+(Yldmk(1)-x_true_1(2,1))^2);
dstlmk(2,1)=sqrt(((Xldmk(2)-x_true_1(1,1))^2)+(Yldmk(2)-x_true_1(2,1))^2);
dstlmk(3,1)=sqrt(((Xldmk(3)-x_true_1(1,1))^2)+(Yldmk(3)-x_true_1(2,1))^2);
    
% true trajectory
for i = 2:N
    % x_k+1 = x_k + V_k * dt * cos(phi_k)
    x_true(1,i) = x_true(1,i-1) + dt*v_true(i-1)*cos(x_true(3,i-1));
    
    % y_k+1 = y_k + V_k *dt * sin(phi_k)
    x_true(2,i) = x_true(2,i-1) + dt*v_true(i-1)*sin(x_true(3,i-1));
    
    % phi_k+1 = phi_k + omega_k * dt
    x_true(3,i) = x_true(3,i-1) + dt*w_true(i-1);
    
     %%rest of distances to landmark ( using standard distance equation)
    dstlmk(1,i)=sqrt(((Xldmk(1)-x_true(1,i))^2)+(Yldmk(1)-x_true(2,i))^2);
    dstlmk(2,i)=sqrt(((Xldmk(2)-x_true(1,i))^2)+(Yldmk(2)-x_true(2,i))^2);
    dstlmk(3,i)=sqrt(((Xldmk(3)-x_true(1,i))^2)+(Yldmk(3)-x_true(2,i))^2);
    
end

% Odometry measurements
v_m = v_true + 0.01*v_true .* randn(1,N-1);
w_m = w_true + 0.04*v_true .* randn(1,N-1);

%%exterceptive measurement
%%distance meaurements (adding noise)
%%This will also be returned
distance = dstlmk+0.01*randn(3,length(dstlmk));
