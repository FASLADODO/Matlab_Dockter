function [x_true, v_m, w_m, bearing] = rws_bear(N, dt, x_true_1, v_true, w_true, Xldmk, Yldmk)
%%Rod Dockter
%%Homework 3
%%problem 3
%%csci 5552

%%real world simulation for bearing only
%%Called by the twoD_EKF in order to calculate the bearing to each landmark
%%and then corrupt it by noise and determine the velocity and angular
%%velocity corrupted by noise and return these values.

% preallocate space
x_true = zeros(3,N);

% initialize
x_true(:,1) = x_true_1;

%%%bearing to landmarks
%%using the -phi+tan^-1(ylandmark-yposition/xlandmark-xposition);
%%Note that all the bearings to all 3 landmarks are calculated even though
%%they arent used
bearlmk = zeros(3,N);
bearlmk(1,1)=-x_true_1(3,1)+atan2(abs(Yldmk(1)-x_true_1(2,1)),abs(Xldmk(1)-x_true_1(1,1)));
bearlmk(2,1)=-x_true_1(3,1)+atan2(abs(Yldmk(2)-x_true_1(2,1)),abs(Xldmk(2)-x_true_1(1,1)));
bearlmk(3,1)=-x_true_1(3,1)+atan2(abs(Yldmk(3)-x_true_1(2,1)),abs(Xldmk(3)-x_true_1(1,1)));

% true trajectory
for i = 2:N
    %%Standard motion update equations
    % x_k+1 = x_k + V_k * dt * cos(phi_k)
    x_true(1,i) = x_true(1,i-1) + dt*v_true(i-1)*cos(x_true(3,i-1));
    
    % y_k+1 = y_k + V_k *dt * sin(phi_k)
    x_true(2,i) = x_true(2,i-1) + dt*v_true(i-1)*sin(x_true(3,i-1));
    
    % phi_k+1 = phi_k + omega_k * dt
    x_true(3,i) = x_true(3,i-1) + dt*w_true(i-1);
    
     %%calculating the rest of bearings to all landmarks
    bearlmk(1,i)=-x_true(3,i)+atan2(abs(Yldmk(1)-x_true(2,i)),abs(Xldmk(1)-x_true(1,i)));
    bearlmk(2,i)=-x_true(3,i)+atan2(abs(Yldmk(2)-x_true(2,i)),abs(Xldmk(2)-x_true(1,i)));
    bearlmk(3,i)=-x_true(3,i)+atan2(abs(Yldmk(3)-x_true(2,i)),abs(Xldmk(3)-x_true(1,i)));
    
end

% Odometry measurements
v_m = v_true + 0.01*v_true .* randn(1,N-1);
w_m = w_true + 0.04*v_true .* randn(1,N-1);

%%exterceptive measurements
%%bearing meaurements (adding noise)
bearing = bearlmk+0.01*randn(3,length(bearlmk));
