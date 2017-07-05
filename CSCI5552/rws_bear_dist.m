function [x_true, v_m, w_m, bearing, distance] = rws_bear_dist(N, dt, x_true_1, v_true, w_true, Xldmk, Yldmk)
%%Rod Dockter
%%Homework 3
%%problem 3
%%csci 5552

%%real world simulation for bearing and distance
%%called in twoD_EKF
%%This calculates both the bearing and distance to all 3 landmarks and then
%%corrupts both independently by noise. It calculates these values using
%%the same equations found in rws and rws_bear.
%%It also calculates the measured velocities and angular velocities

% preallocate space
x_true = zeros(3,N);
% initialize
x_true(:,1) = x_true_1;

%%%bearing to landmarks
%%using the -phi+tan^-1(ylandmark-yposition/xlandmark-xposition);
%%(still for all three landmarks)
bearlmk = zeros(3,N);
bearlmk(1,1)=-x_true_1(3,1)+atan2(abs(Yldmk(1)-x_true_1(2,1)),abs(Xldmk(1)-x_true_1(1,1)));
bearlmk(2,1)=-x_true_1(3,1)+atan2(abs(Yldmk(2)-x_true_1(2,1)),abs(Xldmk(2)-x_true_1(1,1)));
bearlmk(3,1)=-x_true_1(3,1)+atan2(abs(Yldmk(3)-x_true_1(2,1)),abs(Xldmk(3)-x_true_1(1,1)));

%%%distance to landmarks
%%distance = sqrt(xlandmark-xposition^2+ylandmark-yposition^2)
%%(still for all three landmarks)
dstlmk = zeros(3,N);
dstlmk(1,1)=sqrt(((Xldmk(1)-x_true_1(1,1))^2)+(Yldmk(1)-x_true_1(2,1))^2);
dstlmk(2,1)=sqrt(((Xldmk(2)-x_true_1(1,1))^2)+(Yldmk(2)-x_true_1(2,1))^2);
dstlmk(3,1)=sqrt(((Xldmk(3)-x_true_1(1,1))^2)+(Yldmk(3)-x_true_1(2,1))^2);

% true trajectory
for i = 2:N
    %%standard motion equations
    % x_k+1 = x_k + V_k * dt * cos(phi_k)
    x_true(1,i) = x_true(1,i-1) + dt*v_true(i-1)*cos(x_true(3,i-1));
    
    % y_k+1 = y_k + V_k *dt * sin(phi_k)
    x_true(2,i) = x_true(2,i-1) + dt*v_true(i-1)*sin(x_true(3,i-1));
    
    % phi_k+1 = phi_k + omega_k * dt
    x_true(3,i) = x_true(3,i-1) + dt*w_true(i-1);
    
    %%rest of bearings to landmarks
    bearlmk(1,i)=-x_true(3,i)+atan2(abs(Yldmk(1)-x_true(2,i)),abs(Xldmk(1)-x_true(1,i)));
    bearlmk(2,i)=-x_true(3,i)+atan2(abs(Yldmk(2)-x_true(2,i)),abs(Xldmk(2)-x_true(1,i)));
    bearlmk(3,i)=-x_true(3,i)+atan2(abs(Yldmk(3)-x_true(2,i)),abs(Xldmk(3)-x_true(1,i)));
    
    %%rest of distances to landmark
    dstlmk(1,i)=sqrt(((Xldmk(1)-x_true(1,i))^2)+(Yldmk(1)-x_true(2,i))^2);
    dstlmk(2,i)=sqrt(((Xldmk(2)-x_true(1,i))^2)+(Yldmk(2)-x_true(2,i))^2);
    dstlmk(3,i)=sqrt(((Xldmk(3)-x_true(1,i))^2)+(Yldmk(3)-x_true(2,i))^2);
    
end

% Odometry measurements
v_m = v_true + 0.01*v_true .* randn(1,N-1);
w_m = w_true + 0.04*v_true .* randn(1,N-1);

%%exterceptive measurements
%%bearing meaurements (adding noise)
bearing = bearlmk+0.01*randn(3,length(bearlmk));

%%distance meaurements (adding noise)
distance = dstlmk+0.01*randn(3,length(dstlmk));
%%both distance and bearing are returned