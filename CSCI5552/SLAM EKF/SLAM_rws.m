function [x_true, v_m, w_m, z] = SLAM_rws(N, dt, x_true_1, v_true, sigma_v, w_true, sigma_w, XL, NL, sigma_d, sigma_theta, dmax)
%Rod Dockter
%CSCI 5552
%Real world simulation

%Skeleton comes from homework 3 solutions
%returns d theta bearing measurements in addition to the

% pre-allocate space
x_true = zeros(3,N);
v_m = zeros(1,N);
w_m = zeros(1,N);

% initialize
x_true(:,1) = x_true_1;

% true trajectory
for i = 2:N
    
    x_true(1,i) = x_true(1,i-1) + v_true(i-1)*dt*cos(x_true(3,i-1));
    x_true(2,i) = x_true(2,i-1) + v_true(i-1)*dt*sin(x_true(3,i-1));
    x_true(3,i) = x_true(3,i-1) + w_true(i-1)*dt;

    for k = 1:NL
        %%getting error on d and theta measurement
        nd = sigma_d*randn(1);
        ntheta = sigma_theta*randn(1);
        %%calculating the d and theta to landmarks
        d = sqrt((XL(2*(k-1)+1) - x_true(1,i))^2 + (XL(2*(k-1)+2) - x_true(2,i))^2);
        theta = -x_true(3,i) + atan2((XL(2*(k-1)+2) - x_true(2,i)),(XL(2*(k-1)+1) - x_true(1,i)));
        %%if certain landmarks are too far away, set them to zero
        if (d > dmax)
            z(2*(k-1)+1,i) = 0;
        else
            z(2*(k-1)+1,i) = d + nd;
        end
        %%applying error
        z(2*(k-1)+2,i) = theta + ntheta;
        
    end

end

% Odometry measurements
v_m = v_true + sigma_v * v_true.*randn(1,N-1);
w_m = w_true + sigma_w * v_true.*randn(1,N-1);

end



