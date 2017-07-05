%%Rod Dockter
%%Homework 3
%%problem 3
%%csci 5552
% 2-D EKF-Localization 
%%Code uses skeleton code for homework 2 solutions with my EKF update
%%functionality added

clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration settings (update/landmarks)

%%set update mode
%%choose between 'distanceupdate', 'bearingupdate' or 'beardistupdate'
updatemode = 'distanceupdate';
%%set number of landmarks;
%%This can be changed from 1..3 without any alteration (extra credit)
nl = 1;

N = 150; % number of timesteps
dt = 1; % sampling time
%%noise factors
sigma_v = 0.01;
sigma_w = 0.04;
sigma_d = 0.01;
sigma_theta = 0.01;
%%Coordinates of landmarks (random values (centered in circle path)
%%although I always make 3, only the number of landmarks slected gets used
%%This is just to make coding easier
Xldmk = -3+randn(1,3);
Yldmk = -3+randn(1,3);

% Trajectory
v_true = 0.2 + .01*randn(1,N-1); % true velocity. Changing this profile allows changing the trajectory
w_true = -.04*ones(1,N-1); % true rotational velocity

% Updates
Msmt_avail = true; % in this Pb., there are no measurements

% Initial conditions
x_true_1 = [0 0 0]'; % initial starting pose (x y phi)
x_hat_1  = [0 0 0]'; % initial estimate
P_1 = zeros(3); % initial covariance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Real World Simulation
% provides x_true, v_m, z_g

%%based on the updatemode selected the code uses a different real world
%%simulation to get either distance, bearing or distance & bearing
%%measurements along with the measured velocity and measured angular
%%velocity
switch updatemode
    case 'distanceupdate'
        [x_true, v_m, w_m, distance] = rws(N, dt, x_true_1, v_true, w_true, Xldmk, Yldmk);
    case 'bearingupdate'
        [x_true, v_m, w_m, bearing] = rws_bear(N, dt, x_true_1, v_true, w_true, Xldmk, Yldmk);
    case 'beardistupdate'
        [x_true, v_m, w_m, bearing, distance] = rws_bear_dist(N, dt, x_true_1, v_true, w_true, Xldmk, Yldmk);
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bookkeeping

% pre-allocate space for certain values we want to keep track of
x_hat_min = zeros(3,N); % state estimate after Propagation
x_hat_plus = zeros(3,N); % state estimate after update
P_min = zeros(3,3,N); % covariance after Propagation
P_plus = zeros(3,3,N); % covariance after update

%%based on the selected update mode and the number of landmarks used the size of the residual and S that
%%needs to be preallocated will vary for distance and bearing seperate,
%%they need to both only have numberoflandmarks rows, but for the ocmbined
%%the res neeeds numberoflandmarks*2 and S needs #LM*2x#LM*2
switch updatemode
    case 'distanceupdate'
        res = zeros(nl,N); % measurement residual
        S = zeros(nl,nl,N); % residual covariance
    case 'bearingupdate'
        res = zeros(nl,N); % measurement residual
        S = zeros(nl,nl,N); % residual covariance
    case 'beardistupdate'
        res = zeros(2*nl,N); % measurement residual
        S = zeros(2*nl,2*nl,N); % residual covariance
end


% initialize those with the right values where appropriate
%%Guesses
x_hat_plus(:,1) = x_hat_1;
P_plus(:,:,1) = P_1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EKF
%propogation is always the same, update depends on mode being used
for i = 2:N
    
   % Propagation
   [x_hat_min(:,i), P_min(:,:,i)] = EKF_propagate(x_hat_plus(:,i-1), P_plus(:,:,i-1), v_m(i-1), w_m(i-1), v_m(i-1)^2*diag([0.01^2 0.04^2]), dt);
   
   % Update
   %%# of landmarks is passed into each update so only the amount chosen
   %%are considered.
   if Msmt_avail % for this problem, we don't have measurements available
        switch updatemode
            %%calling the various different EKF updates
            case 'distanceupdate'
                [x_hat_plus(:,i), P_plus(:,:,i), res(1:nl,i), S(1:nl,1:nl,i)] = EKF_update_dist(x_hat_min(:,i), P_min(:,:,i), Xldmk, Yldmk, distance(1:nl,i), nl);
            case 'bearingupdate'
                [x_hat_plus(:,i), P_plus(:,:,i), res(1:nl,i), S(1:nl,1:nl,i)] = EKF_update_bear(x_hat_min(:,i), P_min(:,:,i), Xldmk, Yldmk, bearing(1:nl,i), nl);
            case 'beardistupdate'
                [x_hat_plus(:,i), P_plus(:,:,i), res(1:2*nl,i), S(1:2*nl,1:2*nl,i)] = EKF_update_bear_dist(x_hat_min(:,i), P_min(:,:,i), Xldmk, Yldmk, bearing(1:nl,i), distance(1:nl,i), nl);
        end 
   else
        % no update, this wont actually happen for this situation
        x_hat_plus(:,i) = x_hat_min(:,i);
        P_plus(:,:,i) = P_min(:,:,i);
   end
    %%End of EKF once loop has run -1 times
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualization

% labels
state = {'x (m)','y (m)','\phi (rad)'};
stateerr = {'x-x_{hat} (m)','y-y_{hat} (m)','\phi-\phi_{hat} (rad)'};
% generate time 
t = 0:dt:(N-1)*dt;

% Plot of True state vs. posterior est
%%With error ellipses using Esha's code
figure('Name','2D Trajectory'); hold on
plot(x_true(1,:), x_true(2,:), 'b')
plot(x_hat_plus(1,:), x_hat_plus(2,:), 'r')
plot(Xldmk(1:nl), Yldmk(1:nl), 'p','MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',10)
for j = 10:10:150
plot_error_ellipse(x_hat_plus(1:2,j),P_plus(1:2,1:2,j))
plot(x_true(1,j), x_true(2,j), 'xb')
plot(x_hat_plus(1,j), x_hat_plus(2,j), 'xr')
end
xlabel('x (m)')
ylabel('y (m)')
legend('True State','Estimate','Landmarks','error ellipse','true point', 'estimate point')
hold off;


% Plot of Error bounded by 3*sigma
figure('Name','Pose Error'); hold on
for i = 1:3
    subplot(3,1,i); plot(t, x_true(i,:)-x_hat_plus(i,:)); hold on;  
    plot(t, 3*sqrt(squeeze(P_plus(i,i,:))),'r'); hold on;
    plot(t,-3*sqrt(squeeze(P_plus(i,i,:))),'r'); hold on;
    legend('State Error','3\sigma - bound')
    xlabel('time (s)')
    ylabel(stateerr{i})
end
  

%%Plotting the residuals for the diffrent update modes. I print them
%%differently for each mode since the size of res will be different
%%landmarks are displayed in green stars
%%the corresponding true and estimate positions for each error ellipse are
%%also printed on the trajectory
switch updatemode
    case 'distanceupdate'
        % % Plot of the Residual from the update
        [rowres,colres]=size(res);
        figure('Name','Residual'); hold on
        for i = 1:rowres
            subplot(rowres,1,i); plot(t, res(i,:)); hold on;
            subplot(rowres,1,i); plot(t, 3*sqrt(squeeze(S(i,i,:))),'r'); hold on;
            subplot(rowres,1,i); plot(t,-3*sqrt(squeeze(S(i,i,:))),'r'); hold on;
            ylabel(['z_' int2str(i) ' - z_{' int2str(i) '_{hat}} (m)'])
        end
        xlabel('time (s)');
        legend('Measurement Residual','3\sigma - bound')
    case 'bearingupdate'
        % % Plot of the Residual from the update
        [rowres,colres]=size(res);
        figure('Name','Residual'); hold on
        for i = 1:rowres
            subplot(rowres,1,i); plot(t, res(i,:)); hold on;
            subplot(rowres,1,i); plot(t, 3*sqrt(squeeze(S(i,i,:))),'r'); hold on;
            subplot(rowres,1,i); plot(t,-3*sqrt(squeeze(S(i,i,:))),'r'); hold on;
            ylabel(['z_' int2str(i) ' - z_{' int2str(i) '_{hat}} (m)'])
        end
        xlabel('time (s)');
        legend('Measurement Residual','3\sigma - bound')
    case 'beardistupdate'
        % % Plot of the Residual from the update
        [rowres,colres]=size(res);
        figure('Name','Residual'); hold on
        for i = 1:rowres
            subplot(rowres/2,2,i); plot(t, res(i,:)); hold on;
            subplot(rowres/2,2,i); plot(t, 3*sqrt(squeeze(S(i,i,:))),'r'); hold on;
            subplot(rowres/2,2,i); plot(t,-3*sqrt(squeeze(S(i,i,:))),'r'); hold on;
            ylabel(['z_' int2str(i) ' - z_{' int2str(i) '_{hat}} (m)'])
        end
        xlabel('time (s)');
        legend('Measurement Residual','3\sigma - bound')
end

