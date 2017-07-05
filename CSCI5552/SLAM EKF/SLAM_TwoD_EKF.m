% Robotics CSci 5552
% 2-D EKF- SLAM
% Rod Dockter
% Homework 4, problem 2

%Code ported from c++ code for project

clear all;
close all;
clc;

% set time steps and sample time
N = 1000; % number of timesteps
dt = 1; % sampling time

% Choose the number of landmarks <nL>
NL = 24;

% Define the system noise
sigma_v = 0.01; % odometry noise multiplier 
sigma_w = 0.01; % odometry noise multiplier

% Define the measurement noise
sigma_d = 0.01; % distance noise 
sigma_theta = 0.01; % bearing noise

% Define the trajectory
v_true = 0.1*ones(1,N-1); % linear velocity
w_true = 0.01*ones(1,N-1); % rotational velocity

% Define the maximum distance sensing range
dmax = 0.1*v_true(1)/w_true(1);

% Initial conditions
x_true_1 = [0 0 0]'; % initial starting point
x_hat_1 = [0 0 0]'; % initial estimate
P_1 = zeros(3,3); % initial covariance

% coordinates of our known landmarks
lm_center = [0 v_true(1)/w_true(1)]'; %center of circle for robot
% using random offsets get circle path data
for k = 1:NL   
    lm_rad = v_true(1)/w_true(1)+rand(1); % LMs are going to be placed inside circle of .4*radius 
    lm_phi = 2*pi*rand(1);
    % and now the landmark coordinates
    XL(((2*(k-1)+1):(2*(k-1)+2)),1) = [lm_center(1) + lm_rad*cos(lm_phi);...
                                       lm_center(2) + lm_rad*sin(lm_phi)];
end

%%Calling rws to get measurements and odometry readings
[x_true, v_m, w_m, z] = SLAM_rws(N, dt, x_true_1, v_true, sigma_v, w_true, sigma_w, XL, NL, sigma_d, sigma_theta, dmax);


% pre-allocate for EKF Variables
x_hat_m = zeros(3,N); % state estimate after Propagation
x_hat_p = zeros(3,N); % state estimate after update
P_m = zeros(3,3,N); % covariance after Propagation
P_p = zeros(3,3,N); % covariance after update

% measurement noise covariance
R = [sigma_d^2,0;...
     0, sigma_theta^2];

% Preallocate posterior state estimate and covariance
x_hat_p(:,1) = x_hat_1;
P_p(:,:,1) = P_1;
ri =0;

Nc = 0;
for i = 2:N
   % propogate the state
   [x_hat_m(:,i), P_m(:,:,i)] = SLAM_EKF_Propogate(x_hat_p(:,i-1), P_p(:,:,i-1), v_m(1,i-1), sigma_v, w_m(1,i-1), sigma_w, dt, Nc);
   
   %dtermine is landmarks measurements are legitimate (in range)
   clear zValid
   nValid(i) = 0;
   for k = 1:NL
       if z(2*(k-1)+1,i) ~= 0
            nValid(i) = nValid(i)+1;
            zValid((2*(nValid(i)-1)+1):(2*(nValid(i)-1)+2),1) = z((2*(k-1)+1):(2*(k-1)+2),i);
       end
   end
   %if our measurements are useful, update using them
   if (nValid(i) > 0)
        [x_hat_p2,P_p2,res2,S2,Nc] = SLAM_EKF_Update(x_hat_m(:,i), P_m(:,:,i), zValid, Nc, R);
        
        %Resizing the state given new measurements
        if size(x_hat_p2,1) > size(x_hat_m,1)
            oldsize = size(x_hat_m,1);
            newsize = size(x_hat_p2,1);
            %resize the covariance matrix too
            P_pold = P_p;
            P_p = zeros(newsize,newsize,N);
            P_p(1:oldsize,1:oldsize,:) = P_pold;
            P_p(:,:,i) = P_p2;
            
            P_mold = P_m;
            P_m = zeros(newsize,newsize,N);
            P_m(1:oldsize,1:oldsize,:) = P_mold;

            x_hat_pold = x_hat_p;
            x_hat_p = zeros(newsize,N);
            x_hat_p(1:oldsize,:) = x_hat_pold;
            x_hat_p(:,i) = x_hat_p2;

            x_hat_mold = x_hat_m;
            x_hat_m = zeros(newsize,N);
            x_hat_m(1:oldsize,:) = x_hat_mold;
        else
            x_hat_p(:,i) = x_hat_p2;
            P_p(:,:,i) = P_p2;
        end
        
        for j = 1:(length(res2)/2)
            res(1:2,ri+1) = res2((2*j-1):(2*j),1);
            resCov(1,ri+1) = S2((2*j-1),(2*j-1));
            resCov(2,ri+1) = S2((2*j),(2*j));
            ri = ri + 1;
        end
        
   else
        x_hat_p(:,i) = x_hat_m(:,i);
        P_p(:,:,i) = P_m(:,:,i);
   end
   
   
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plotting%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%stuff%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%Plotting each individual trajectory vs actual
% generate time
t = 0:dt:(N-1)*dt;

% True state vs. posterior est
figure('Name','Individual State Trajectory');
subplot(3,1,1)
    plot(t, x_true(1,:), 'b',t, x_hat_p(1,:), 'r');
    xlabel('Time (s)')
    ylabel('X Position (m)')
    legend('True State','Estimated State','Location','NorthWest')
subplot(3,1,2)
    plot(t, x_true(2,:), 'b',t, x_hat_p(2,:), 'r');
    xlabel('Time (s)')
    ylabel('Y Position (m)')
    legend('True State','Estimated State','Location','NorthWest')
subplot(3,1,3)
    plot(t, x_true(1,:), 'b',t, x_hat_p(1,:), 'r');
    xlabel('Time (s)')
    ylabel('Orientation (radrees)')
    legend('True State','Estimated State','Location','NorthWest')
%%
% Plotting the total 2-d trajectory with the acutal trajectory and the
% landmarks

% using "plot_error_ellipse.m"
npts = 100; 
th = linspace(0, 2*pi, npts);
r = 3;

figure('Name','Trajectory');
plot(x_true(1,:),x_true(2,:),'b',x_hat_p(1,:),x_hat_p(2,:),'r');
    title('Trajectory') 
    xlabel('X Position (m)')
    ylabel('Y Position (m)')
    hold on
    for i = 1:NL
        plot(XL(2*(i-1)+1),XL(2*(i-1)+2),'h',...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10)
    end
    for i = 1:Nc
        plot(x_hat_p(2*(i-1)+4,end),x_hat_p(2*(i-1)+5,end),'h',...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','b',...
            'MarkerSize',5)
        
        [Evec, Eval] = eig(P_p((2*(i-1)+4):(2*(i-1)+5),(2*(i-1)+4):(2*(i-1)+5),end));
        ptsl = repmat(x_hat_p((2*(i-1)+4):(2*(i-1)+5),end),1,npts) + r*Evec*sqrt(Eval)*[cos(th); sin(th)];
        plot(ptsl(1,:), ptsl(2,:),'--k');
    end
    for i = 1:20
         % Used "plot_error_ellipse.m"
        [Evec, Eval] = eig(P_p(1:2,1:2,(i-1)*N/20+1));
        pts = repmat(x_hat_p(1:2,(i-1)*N/20+1),1,npts) + r*Evec*sqrt(Eval)*[cos(th); sin(th)];
        plot(x_true(1,(i-1)*N/20+1),x_true(2,(i-1)*N/20+1),'xb');
        plot(x_hat_p(1,(i-1)*N/20+1),x_hat_p(2,(i-1)*N/20+1),'xr');
        plot(pts(1,:), pts(2,:),'--k');
    end

    hold off
    legend('True State','Estimated State','Landmarks','Location','NorthWest')
%%
% Plotting the error of 3 different poses with the 3 sigma bounds
figure('Name','Estimate Errors');
subplot(3,1,1)
    plot(t, x_true(1,:)-x_hat_p(1,:), 'b',...
        t, 3*sqrt(shiftdim(P_p(1,1,:))),'r',...
        t,-3*sqrt(shiftdim(P_p(1,1,:))),'r');
    title('X Position Estimate Error')
    xlabel('Time (s)')
    ylabel('Error (m)')
    legend('State Estimate Error','3\sigma - bound','Location','NorthWest')
subplot(3,1,2)
    plot(t, x_true(2,:)-x_hat_p(2,:), 'b',...
        t, 3*sqrt(shiftdim(P_p(2,2,:))),'r',...
        t,-3*sqrt(shiftdim(P_p(2,2,:))),'r');
    title('Y Position Estimate Error')    
    xlabel('Time (s)')
    ylabel('Error (m)')
    legend('State Estimate Error','3\sigma - bound','Location','NorthWest')
subplot(3,1,3)
    plot(t, x_true(3,:)-x_hat_p(3,:), 'b',...
        t, 3*sqrt(shiftdim(P_p(3,3,:))),'r',...
        t,-3*sqrt(shiftdim(P_p(3,3,:))),'r');
    title('Orientation Estimate Error')    
    xlabel('Time (s)')
    ylabel('Error (rad)')
    legend('State Estimate Error','3\sigma - bound','Location','NorthWest')
%%
% Error
% Error
figure('Name','Estimate Errors');
subplot(3,1,1)
    plot(t, x_true(1,:)-x_hat_p(1,:), 'b',...
        t, 3*sqrt(shiftdim(P_p(1,1,:))),'r',...
        t,-3*sqrt(shiftdim(P_p(1,1,:))),'r');
    title('X Position Estimate Error')
    xlabel('Time (s)')
    ylabel('Error (m)')
    legend('State Estimate Error','3\sigma - bound','Location','NorthWest')
subplot(3,1,2)
    plot(t, x_true(2,:)-x_hat_p(2,:), 'b',...
        t, 3*sqrt(shiftdim(P_p(2,2,:))),'r',...
        t,-3*sqrt(shiftdim(P_p(2,2,:))),'r');
    title('Y Position Estimate Error')    
    xlabel('Time (s)')
    ylabel('Error (m)')
    legend('State Estimate Error','3\sigma - bound','Location','NorthWest')
subplot(3,1,3)
    plot(t, x_true(3,:)-x_hat_p(3,:), 'b',...
        t, 3*sqrt(shiftdim(P_p(3,3,:))),'r',...
        t,-3*sqrt(shiftdim(P_p(3,3,:))),'r');
    title('Orientation Estimate Error')    
    xlabel('Time (s)')
    ylabel('Error (rad)')
    legend('State Estimate Error','3\sigma - bound','Location','NorthWest')

% Residuals
figure('Name','Residuals');
subplot(2,1,1)
    plot(1:ri,res(1,:), 'b',...
        1:ri, 3*sqrt(resCov(1,:)),'r',...
        1:ri,-3*sqrt(resCov(1,:)),'r');
    title('X Position Measurement Estimate Error')
    xlabel('Measurement')
    ylabel('Error (m)')
    legend('Residual','3\sigma - bound','Location','NorthWest')
subplot(2,1,2)
    plot(1:ri,res(2,:), 'b',...
        1:ri, 3*sqrt(resCov(2,:)),'r',...
        1:ri,-3*sqrt(resCov(2,:)),'r');
    title('Y Position Measurement Estimate Error')
    xlabel('Measurement')
    ylabel('Error (m)')
    legend('Residual','3\sigma - bound','Location','NorthWest')