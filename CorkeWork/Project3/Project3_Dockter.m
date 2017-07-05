%% Project 3
% Rod Dockter
% ME 5286


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part A - Free Motion

%pos,vel,acc
q0 = [0 +pi/2 -pi/2 0 0 0];     
qd0 = -[1 1 1 1 1 1];   % initial robot joint velocities (invoked within
                        % s1_ztorqueFreefalls's Robot blck
                        % (within it's 1st internal integrator)
q_in = [ 0 0 0 0 0 0]'; % Torque input to motors; this is what you'll want 
                        % to change in a computed torque controller

% SeUsing Puma
mdl_puma560;
myRobot = p560;  % myRobot is used internally by sl_ztroqueFreefall
myRobot = p560.nofriction;  % drop coloumb friction for faster sim times
                            % This still includes viscous friction



% simulate the robot in free fall
figure(1); clf
open('sl_ztorqueFreefall.slx');
r = sim('sl_ztorqueFreefall', 'StopTime', '10' )  ; % run simulation

% get info out of simulink
t = r.get('tout')';
q_total = r.get('yout')';
%pos,vel,acc
q = q_total(1:6,:);
qdot = q_total(7:12,:);
qdotdot = q_total(13:18,:);


[row,col] = size(q);

for(kk = 1:col)
    tempP = myRobot.fkine(q(1:6,kk)');
    CartPos(kk) = norm(tempP(1:3,4));
    J = myRobot.jacob0(q(1:6,kk)');
    tempV = J*qdot(1:6,kk);
    CartVel(kk) = norm(tempV);
    tempA = J*qdotdot(1:6,kk);
    CartAcc(kk) = norm(tempA);
end

% plot stuff
figure(2); clf
plot(t,q);
xlabel('Time [s]')
ylabel('Joint Angles [rad]')
legend([['qqqqqq']' ['123456']' ] )
title('\bfRobot Free Fall')

figure(3); clf
plot(t,CartPos,'r');
hold on
plot(t,CartVel,'b');
hold on
plot(t,CartAcc,'g');
hold off
xlabel('Time [s]')
ylabel('Cartesian Norms')
legend('Position','Velocity','Acceleration' )
title('\bfRobot Free Fall')






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part B
% Let's implement gravity compensation (LETS NOT)

%straight up
q0 = [0 +pi/2 -pi/2 0 0 0];     
qd0 = -[1 1 1 1 1 1];   % initial robot joint velocities (invoked within
                        % s1_ztorqueFreefalls's Robot blck
                        % (within it's 1st internal integrator)
q_in = [ 0 0 0 0 0 0]'; % Torque input to motors; this is what you'll want 
                        % to change in a computed torque controller
 

TauG = p560.gravload(q0)


qd  = 0*q0;  % ignoring velocity effects
qdd = 0*q0;  % ignoring acceleration effects
Tau = myRobot.rne(q0, 0*q0, 0*q0)   % , GRAV, FEXT) extras:


% simulate the robot in free fall WITH gravity compensation
% we'll use simulink anyway just to get familiar with it.
figure(1); clf; % Corke's toolkit can only plot single robot, use same fig
open('sl_ztorqueGravComp.slx');
r = sim('sl_ztorqueGravComp', 'StopTime', '10' )  ;       

% get info out of simulink
t_g = r.get('tout')';
q_g = r.get('yout')';

q = q_g(1:6,:);
qdot = q_g(7:12,:);
qdotdot = q_g(13:18,:);

[row,col] = size(q_g);

for(kk = 1:col)
    tempP = myRobot.fkine(q(1:6,kk)');
    CartPos(kk) = norm(tempP(1:3,4));
    J = myRobot.jacob0(q(1:6,kk)');
    tempV = J*qdot(1:6,kk);
    CartVel(kk) = norm(tempV);
    tempA = J*qdotdot(1:6,kk);
    CartAcc(kk) = norm(tempA);
end

% plot stuff
figure(2); clf
plot(t_g,q);
xlabel('Time [s]')
ylabel('Joint Angles [rad]')
legend([['qqqqqq']' ['123456']' ] )
title('\bfRobot Gravity Compensation')

figure(3); clf
plot(t_g,CartPos,'r');
hold on
plot(t_g,CartVel,'b');
hold on
plot(t_g,CartAcc,'g');
hold off
xlabel('Time [s]')
ylabel('Cartesian Norms')
legend('Position','Velocity','Acceleration' )
title('\bfRobot Gravity Compensation')

%-----Given viscous friction and no commanded inputs: 
%-----Is a gravity-compensated robot stable to impulse disturbances at the end effector?
%Yes by inspection the joints fall initally and then stabilize

% -----Would the gravity compensation commanded torques change if the 
% -----Mass/Inertia matrix had significant errors?
% Yes because T = I*alpha duh
% Kinematics had significant errors?
% Yes because the mass matrix depends on the FK

%% Part B cont.

%-----Provide all joint configurations for when the robot joint torques are un-affected by gravity
%-----(show a 3D plot or sketch of the robot’s pose similar to this): 
%There are four positions by inspection:

q_stab_1 = [0 +pi/2 -pi/2 0 0 0]; 
q_stab_2 = [0 +pi/2 +pi/2 0 0 0]; 
q_stab_3 = [0 -pi/2 -pi/2 0 0 0]; 
q_stab_4 = [0 -pi/2 +pi/2 0 0 0]; 

figure(1)
myRobot.plot(q_stab_1);
title('straight up')

figure(2)
myRobot.plot(q_stab_2);
title('straight up, elbow down')

figure(3)
myRobot.plot(q_stab_3);
title('straight down')

figure(4)
myRobot.plot(q_stab_4);
title('straight down, elbow up')

%% Wus dis stuff
% figure(5); clf
% p   = myRobot.fkine(q'        );     % get FK
% p_g = myRobot.fkine(q_g'      );
% p   = squeeze(p(1:3, end, :)  );     % extract positions only
% p_g = squeeze(p_g(1:3, end, :));     % squeeze drops unused dimensions
% 
% plot3( p(1,:)  , p(2,:)  , p(3,:)  , 'r.-'); hold on;
% plot3( p_g(1,:), p_g(2,:), p_g(3,:), 'b.-'); grid on
% legend('no gravity compensation', 'gravity compensation')
% xlabel('x [cm]')
% ylabel('y [cm]')
% zlabel('z [cm]')
% title('\bfTask Space End Effector Motion')






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part C - Trajectory Simulate

%TIMS SIMULINK IS BAD!
%Had to adjust to actually do IK

figure(1); clf;
open('sl_ctorqueTraj.slx');  % Different from Tims model (actually does IK)

r = sim('sl_ctorqueTraj', 'StopTime', '2' )  ;       

% get info out of simulink
t_out = r.get('tout')';
y = r.get('yout')';
q = y(     [1:6] ,  :   );
e = y(     [7:12] ,  :   );
%Get XYZ and qd from output nodes in simulink
x_cmd = y(     [13] ,  :   );
y_cmd = y(     [14] ,  :   );
z_cmd = y(     [15] ,  :   );
qd_cmd = y(     [16:21] ,  :   );
qd_act = y(     [22:27] ,  :   );

% plot stuff
figure(40); clf
plot(t_out,q, '.'); hold on
plot(t_out,e,'--');
xlabel('Time [s]')
ylabel('Joint Angles [rad]')
legend([['qqqqqq' 'eeeeee']' ['123456' '123456']' ] )
title('\bfRobot Trajectory Control')

% Joint errors go to zero right  away

%% Part C plot trajectorys
figure(46); clf

% actual 
myTraj2 = myRobot.fkine( q' );
[X,Y,Z] = transl(myTraj2);
plot3( X, Y, Z, 'b.-'); hold on; grid on

plot3( x_cmd, y_cmd, z_cmd, 'r.-'); hold on; grid on

xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
axis equal % maintain aspect ratio


legend ('Actual Path', 'Desired Path')

% They appear right on top of each other!

%% Part C Answers

% ------ Explain your rationale for where you will place the circle relative to the base
%- I chose a 0.05 radius, circle rotated by x and y but with no offsets
%- This gave good trajectory motion without hitting a singularity

% ------ Explain your rationale for controller decisions, design, and gains.
% - I chose to use workspace control so that it was easier to specify my
% trajectory.
% - Then I used torque based controls on the joints, this seemed more ideal
% from reading Corke 9.4.3.2
% - For gains I used Kp = 500 and Kd = 100. To find these I ramped up the
% gains until my simulation seemed to run smoothly.


% ----- Show a plot of the magnitude of error vs time for the entire trajectory, indicating that you meet 
% the position spec and final trajectory time.  Plot error for both position and velocity in the task
% space.

for kk = 1:length(t_out)
    error(kk) = sqrt((X(kk) - x_cmd(kk))^2 + (Y(kk) - y_cmd(kk))^2 + (Z(kk) - z_cmd(kk))^2 );
end


figure(108)
plot( t_out, error, 'r.-');
title('Position error vs. time')
xlabel('t [s]')
ylabel('error [m]')

% Position error settles to below 0.5 mm after an initial spike


% ----- Plot of velocity error vs time

for(kk = 1:length(t_out))
    J = myRobot.jacob0(q(:,kk)');
    tempV = J*qd_cmd(:,kk);
    Vel_Cmd(kk) = norm(tempV(1:3));
    tempV2 = J*qd_act(:,kk);
    Vel_act(kk) = norm(tempV2(1:3));
end


for kk = 1:length(t_out)-1
    errorv(kk) = sqrt( (Vel_Cmd(kk) - Vel_act(kk))^2);
end

figure(109)
plot( t_out(1:end-1), errorv, 'r.-');
title('Velocity error vs. time')
xlabel('t [s]')
ylabel('error [m/s]')

% Velocity error  is around 0.001 m/s after and initial spike which is to
% be expected since the robot starts from rest

