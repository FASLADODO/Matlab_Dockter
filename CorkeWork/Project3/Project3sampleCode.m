%% Project 3: Dynamic Control



%% Constants, mostly for first two simulink files
q0 = [0 0 0 0 0 0];     % initial robot pose (invoked within 
                        % sl_ztroqueFreefall's Robot block 
                        % (within it's second internal integrator)
qd0 = -[1 1 1 1 1 1];   % initial robot joint velocities (invoked within
                        % s1_ztorqueFreefalls's Robot blck
                        % (within it's 1st internal integrator)
q_in = [ 0 0 0 0 0 0]'; % Torque input to motors; this is what you'll want 
                        % to change in a computed torque controller

% Select a model
mdl_puma560;
myRobot = p560;  % myRobot is used internally by sl_ztroqueFreefall
myRobot = p560.nofriction;  % drop coloumb friction for faster sim times
                            % This still includes viscous friction

% % Select a different model
% mdl_stanford;
% myRobot = stanf.nofriction;  % myRobot is used internally by sl_ztroqueFreefall


% simulate the robot in free fall
figure(1); clf
open('sl_ztorqueFreefall.slx');
r = sim('sl_ztorqueFreefall', 'StopTime', '5' )  ; % run simulation

% get info out of simulink
t = r.get('tout')';
q = r.get('yout')';

% plot stuff
figure(2); clf
plot(t,q);
xlabel('Time [s]')
ylabel('Joint Angles [rad]')
legend([['qqqqqq']' ['123456']' ] )
title('\bfRobot Free Fall')


%% Let's implement gravity compensation 
% what does each joint torque have to be to counteract gravity at q = q0?
TauG = p560.gravload(q0)

% look into this, i.e., type edit gravload and review the actual code,
% it's really just RNE with zero joint vel. and acceleration:
qd  = 0*q0;  % ignoring velocity effects
qdd = 0*q0;  % ignoring acceleration effects
Tau = myRobot.rne(q0, 0*q0, 0*q0)   % , GRAV, FEXT) % extras: Grav redefines 
                                    % gravity vector, Fext defines forces 
                                    % applied to End Effector
% The lesson here:  You don't need simulink for this problem.


% simulate the robot in free fall WITH gravity compensation
% we'll use simulink anyway just to get familiar with it.
figure(1); clf; % Corke's toolkit can only plot single robot, use same fig
open('sl_ztorqueGravComp.slx');
r = sim('sl_ztorqueGravComp', 'StopTime', '5' )  ;       

% get info out of simulink
t_g = r.get('tout')';
q_g = r.get('yout')';

% plot stuff
figure(4); clf
plot(t_g,q_g);
xlabel('Time [s]')
ylabel('Joint Angles [rad]')
legend([['qqqqqq']' ['123456']' ] )
title('\bfRobot with Gravity Compensation')

%%
figure(5); clf
p   = myRobot.fkine(q'        );     % get FK
p_g = myRobot.fkine(q_g'      );
p   = squeeze(p(1:3, end, :)  );     % extract positions only
p_g = squeeze(p_g(1:3, end, :));     % squeeze drops unused dimensions

plot3( p(1,:)  , p(2,:)  , p(3,:)  , 'r.-'); hold on;
plot3( p_g(1,:), p_g(2,:), p_g(3,:), 'b.-'); grid on
legend('no gravity compensation', 'gravity compensation')
xlabel('x [cm]')
ylabel('y [cm]')
zlabel('z [cm]')
title('\bfTask Space End Effector Motion')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here is the Cartesian trajectory you want to control:
syms t A k P Pd Pdd

x =  A*sin(k*t);
y =  A*cos(k*t); 
z = 0; 

% this symbolic form might come in handy later
P   = [  x y z 0 ]   *  trotx( -pi/4 )  * troty( pi/4 );
Pd  = diff(  P , t);
Pdd = diff(  Pd, t);

% substitute in real numbers
t = [0:200]/100;   % time in sec
A = 0.05;         % circle radius in meters
k = 2*pi;         % angular freq. < ---- This changes your speed.

myTraj = zeros( length(t), 3);
for j = 1:3
    myTraj(:,j) = eval(P(j));
end

% plot the target directory:
figure(45); clf
plot3( myTraj(:,1), myTraj(:,2), myTraj(:,3), 'r.-'); hold on; grid on

xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title('\bfDesired Trajectory')
axis equal % maintain aspect ratio

%% Now Control the robot (you need to fix this so it does the correct trajectory):
% Corke examples (read the Corke & Sastry sections very closely!)
%sl_fforward; % does Feedforward Control 9.4.3.1 Corke (PD Control Sastry, 4.5.3)
%sl_ctorque;  % does computed torque control 9.4.3.1 Corke (PD Control Sastry, 4.5.3)

%% ME8287 Simulink Models for trajectory following ...
% Corke's examples aren't great for arbitrary trajectories.  
% Prof. K made a trajectory-generating function that is compatible with 
% simulink and modified the simulink model to use it.
% You need to edit that simulink function to invoke your own trajectory

figure(1); clf; % Corke's toolkit can only plot single robot, use same fig
open('sl_ctorqueTraj.slx');  % Different model, uses 'myRobot' and custom
                             % trajectory generation:  desiredTrajectory()
                             % note that it's a simulink function so it'll
                             % compile to mex and you'll need to modify it 
                             % within Simulink to alter it
r = sim('sl_ctorqueTraj', 'StopTime', '5' )  ;       

% get info out of simulink
t = r.get('tout')';
y = r.get('yout')';
q = y(     [1:6] ,  :   );
e = y( 6 + [1:6] ,  :   );

% plot stuff
figure(40); clf
plot(t,q, '.'); hold on
plot(t,e,'--');
xlabel('Time [s]')
ylabel('Joint Angles [rad]')
legend([['qqqqqq' 'eeeeee']' ['123456' '123456']' ] )
title('\bfRobot Trajectory Control')



%% Plot both desired traj and actual traj
figure(46); clf

% actual 
myTraj2 = myRobot.fkine( q' );
[X,Y,Z] = transl(myTraj2);
plot3( X, Y, Z, 'b.-'); hold on; grid on

% desired
plot3( myTraj(:,1), myTraj(:,2), myTraj(:,3), 'r.-'); hold on; grid on
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
axis equal % maintain aspect ratio


legend ('Actual Path', 'Desired Path')







