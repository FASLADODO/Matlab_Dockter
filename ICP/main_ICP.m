% Main Program to test laser scan matching with ICP

close all
clear all
clc

% -------------------------------------------------------------------------
% Parameters
sigma_r = 0.01; % laser distance measurement noise std (m)
sigma_th = 1e-5; % laser bearing msmt. noise std (rad)



% -------------------------------------------------------------------------
% generate two scans from random poses in the environment

% pose 1 (ground truth)
p1 = [2 3]'; % position in global frame (m)
phi1 = pi/4; % bearing in global frame (rad)

% and generate the matching laser scan
% scan1 is a struct with fields scan1.r (true distance), 
% scan1.rm (measured distance), and scan.theta (bearing)
% inside this function, you can specify the obstacles and other params
% inside, set the visualize-flag to 1 to see what's going on...
% notice that if scan1.rm(i)=0 then the laser did not get any valid return
scan1 = generate_scan(p1(1), p1(2), phi1, sigma_r);

% pose 2 (ground truth)
p2 = [2.2 3.1]';
phi2 = pi/4 + 5/180*pi;

scan2 = generate_scan(p2(1), p2(2), phi2, sigma_r);


% -------------------------------------------------------------------------
% Find ground-truth (and initial estimate) for R1_p_R2 and R1_phi_R2

% Relative Position ( from p2 = p1 + C(phi1) * R1_p_R2 )
R1_p_R2 = [ cos(phi1) sin(phi1)
           -sin(phi1) cos(phi1) ] * (p2 - p1);
       
% Relative orientation       
R1_phi_R2 = phi2 - phi1; 
R1_phi_R2 = atan2(sin(R1_phi_R2), cos(R1_phi_R2)); % limit between -pi and pi


% Form initial guess (perturb these values slightly)
R1_p_R2_ini = R1_p_R2 + 0.1*randn(2,1); % (roughly +/- 30 cm 3sigma)
R1_phi_R2_ini = R1_phi_R2 + 0.02*randn; % (roughly +/- 3.5 deg 3sigma)



% -------------------------------------------------------------------------
% Based on the above, find ICP WLS estimate for relative pose
%
% You have at your disposal functions to
% - compute R1_p_li = pol2cart(scan.theta(i), scan.r(i))
% - compute R1_Ri = get_Ri(scan.r(i), scan.theta(i), sigma_r, sigma_theta)
%
% You should implement the following function, which takes the initial
% guess and the two scans, and computes the WLS estimate for the relative
% pose, plus its covariance
% [R1_p_R2_hat, R1_phi_R2_hat, P] = icp_wls(R1_p_R2_ini, R1_phi_R2_ini, scan1, scan2, sigma_r, sigma_th);
%
% You should also implement a function that does ICP matching, i.e.,
% assigns points from scan2 that correspond to the same point measured in
% scan1, based on physical (or Mahalanobis) distance after an estimated
% transform (more detail in HW assignment). This function will likely be
% called by icp_wls().
% [match1, match2] = ICP(R1_p_R2_ini, R1_phi_R2_ini, scan1, scan2, sigma_r, sigma_th);
% where match1 and match2 are vectors of the same length, and
% scan1.r(match1(i)) is the distance to the same point as
% scan2.r(match2(i))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start of work in here

%%Assignment 2, Problem 4
%%This portion of code will act to call ICP_match, pass those parameters into
%%ICP_wls and then run through the iterative looping of the IWLS method.


%%Initializing the xhat estimate
xhat_vector(:,1) = [R1_p_R2_ini; R1_phi_R2_ini];
%%record last iteration step
step =0;

for k = 1:100;
    step = k+1;
    substep = k;
    %%Calling to the weighted least squares subroutine
    [R1_p_R2_hat, R1_phi_R2_hat, P_xx_inv] = ICP_wls(R1_p_R2_ini,R1_phi_R2_ini,scan1,scan2,sigma_r,sigma_th);
    %%building the 3x1 xhat vector
    xhat_vector(:,k+1) = [R1_p_R2_hat; R1_phi_R2_hat];
    %%checking the cost function for a tolerance
    costfunc(1,k) = norm(xhat_vector(:,k+1)-xhat_vector(:,k))/norm(xhat_vector(:,k));
    if costfunc(1,k) < 0.001
       break; 
    end
    %%Updating the initial estimates for the next iteration
    R1_p_R2_ini = R1_p_R2_hat;
    R1_phi_R2_ini = R1_phi_R2_hat;
end

%%inverting to get covariance
P_xx=inv(P_xx_inv);
%%Calculating error
x_error = zeros(3,step);
for j = 1:step
   x_error(:,j) = [R1_p_R2;R1_phi_R2] - xhat_vector(:,j);
end
    
R1_p_R2_hat
R1_phi_R2_hat

%%display position orientation error
figure('Name','Position and Orientation error');
subplot(3,1,1);
    plot(1:step, x_error(2,:), '.b')
    xlabel('Iteration step')
    ylabel('X Position Error')
    legend('xtrue - xhat')
subplot(3,1,2);
    plot(1:step, x_error(2,:), '.b')
    xlabel('Iteration step')
    ylabel('Y Position Error')
    legend('ytrue - yhat')
subplot(3,1,3);
    plot(1:step, x_error(3,:), '.b')
    xlabel('Iteration step')
    ylabel('Phi pose Error')
    legend('phitrue - phihat')
%%Display cost function figure
figure('Name','Cost Function Value');
    plot(1:substep, costfunc(1,:), '.b')
    xlabel('Iteration step')
    ylabel('Cost Function')
    legend('Cost Function')



