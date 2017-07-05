% Robotics CSci 5552
% 1-D EKF-Localization with GPS updates


clear all;
close all;
clc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration settings

N = 150; % number of timesteps
dt = 1; % sampling time

% Trajectory
sigma_v = 0.2; % odometry noise std (not covariance!)
v_true = 0.2 + .01*randn(1,N-1); % true velocity. Changing this profile allows changing the trajectory

% Updates
GPS_avail = ones(1,N); % flag indicating if GPS msmt. is available at timestep i
sigma_g = 0.1; % GPS measurement noise std (not covariance!)

% Initial conditions
x_true_1 = 0; % initial starting point
x_hat_1 = 0; % initial estimate
P_1 = 0; % initial covariance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Real World Simulation
% provides x_true, v_m, z_g

[x_true, v_m, z_g] = rws(N, dt, x_true_1, v_true, sigma_v, sigma_g);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bookkeeping

% pre-allocate space for certain values we want to keep track of
x_hat_min = zeros(1,N); % state estimate after Propagation
x_hat_plus = zeros(1,N); % state estimate after update
x_hat_odo = zeros(1,N); % state estimate without updates (dead reckoning)
P_odo = zeros(1,N); % covariance of dead reckoning only
P_min = zeros(1,N); % covariance after Propagation
P_plus = zeros(1,N); % covariance after update
res = zeros(1,N); % measurement residual
S = zeros(1,N); % residual covariance

% initialize those with the right values where appropriate
x_hat_plus(1,1) = x_hat_1;
x_hat_odo(1,1) = x_hat_1;
P_plus(1,1) = P_1;
P_odo(1,1) = P_1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EKF

% notice that we let the filter start with a propagation step. The first
% GPS measurement at timestep i=1 is discarded (one could implement it
% differently)
for i = 2:N
    
   % Propagation
   [x_hat_min(1,i), P_min(1,i)] = EKF_propagate(x_hat_plus(1,i-1), P_plus(1,i-1), v_m(1,i-1), sigma_v, dt);
   
   % Calculate the dead-reckoning estimate (trajectory based purely on
   % odometry, no update
   [x_hat_odo(1,i), P_odo(1,i)] = EKF_propagate(x_hat_odo(1,i-1), P_odo(1,i-1), v_m(1,i-1), sigma_v, dt);
   
   
   % Update
   if GPS_avail(1,i)
       [x_hat_plus(1,i), P_plus(1,i), res(1,i), S(1,i)] = EKF_update(x_hat_min(1,i), P_min(1,i), z_g(1,i), sigma_g);
   end
    
    
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualization

%%Rod Dockter
%%Assignment 2 problem 2

% generate time and "double time" (for showing prior and posterior
% estimate)
t = 0:dt:(N-1)*dt;
tt = [0 0]; for i = 1:N-1, tt = [tt i*dt*[1 1]]; end;

% generate error that contains both values after propagation and after
% update (to generate saw-tooth-pattern)
ddx = []; for i = 1:N, ddx = [ddx x_true(1,i)-x_hat_min(1,i)  x_true(1,i)-x_hat_plus(1,i)]; end;
ddP = []; for i = 1:N, ddP = [ddP P_min(1,i) P_plus(1,i)]; end;


% True state vs. posterior est
figure('Name','Trajectory'); hold on
plot(t, x_true, 'b'); 
plot(t, x_hat_plus, 'r');
plot(t, x_hat_odo, 'g');
xlabel('time (s)')
ylabel('Position (m)')
legend('True State','Posterior Estimate','Dead Reckoning')


% Error
figure('Name','Dead Reckoning Error'); hold on
plot(t, x_true - x_hat_odo);
plot(t, 3*sqrt(P_odo),'r')
plot(t,-3*sqrt(P_odo),'r')
xlabel('time (s)')
ylabel('x-x_{hat} (m)')
legend('Dead Reckoning Position Error','3\sigma - bound')


% Error
figure('Name','Position Error'); hold on
% plot(t, x_true-x_hat_plus)  % This part will show only posterior error
% plot(t, 3*sqrt(P_plus),'r')
% plot(t,-3*sqrt(P_plus),'r')
plot(tt, ddx) % this part shows saw-tooth pattern
plot(tt, 3*sqrt(ddP),'r')
plot(tt,-3*sqrt(ddP),'r')
xlabel('time (s)')
ylabel('x-x_{hat} (m)')
legend('Position Error with GPS updates','3\sigma - bound')

% Residual
figure('Name','Residual'); hold on
plot(t, res)
plot(t, 3*sqrt(S),'r')
plot(t,-3*sqrt(S),'r')
xlabel('time (s)');
ylabel('z - z_{hat} (m)')
legend('GPS Residual','3\sigma - bound')

