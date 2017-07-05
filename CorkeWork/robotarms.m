%% Robotics arms

clear all

% dh parameters
theta = 0; %no angle offset
d = 0.1; %along z
a = 0.2; % along x
alpha = pi/2; %rotate x
joint = 0; %revolute
L1 = Link([theta,d,a,alpha,joint]);

L1.A(0.5)

%% 2 link arm

clear all

L(1) = Link([0 0 1 0]); %link 1
L(2) = Link([0 0 1 0]); %link 2

two_link = SerialLink(L, 'name', 'two_link');

mdl_twolink

two_link.plot([0,0])

%% gregory

clear all

%offset zero pose
theta_offset = [0,-pi/2,pi/2,0,0,0]; %normally would be in DH table

%theta,d,a,alpha (DH convention)
L(1) = Link([0 0 0 -pi/2]); %link 1
L(2) = Link([0 0 10 0]); %link 2
L(3) = Link([0 0 0 pi/2]); %link 3
L(4) = Link([0 11 0 -pi/2]); %link 4
L(5) = Link([0 0 0 pi/2]); %link 5
L(6) = Link([0 0 0 0]); %link 6 

%create robot
greg = SerialLink(L, 'name', 'gregory');
greg.tool = transl(0,0,9); %end effector offset

% [blank, theta1, theta2, theta3, theta4, theta5,theta6];
qn = [0,0,0,0,0,0] %straight up
qe = [0,0,pi/2,0,0,0]; %elbow out
qw = [0,0,pi/2,0,-pi/2,0]; %elbow out, wrist up
qs = [0,pi/2,-pi/2,0,pi/2,0]; %shoulder over

%Plot particular positions
greg.plot(qn + theta_offset)
% pause(4);
% greg.plot(qe + theta_offset)
% pause(4);
% greg.plot(qw + theta_offset)
% pause(4);
% greg.plot(qs + theta_offset)

% % sweep all angles
% for t=0:.2:16.3
%     greg.plot([1,1,1,1,1,1]*t + theta_offset)
%     pause(1/60);
% end

%move through positions
T0 = greg.fkine(qn+theta_offset)
T1 = greg.fkine(qe+theta_offset);
T2 = greg.fkine(qw+theta_offset);
T3 = greg.fkine(qs+theta_offset);

% closed form solution DOES NOT WORK WELL
qi = greg.ikine6s(T2,'lu')

pause(2);

greg.plot(qi)


% q1 = greg.ikine(T0,theta_offset); % Numerical solution
% q2 = greg.ikine(T1,theta_offset); % Numerical solution
% q3 = greg.ikine(T2,q2); % Numerical solution
% q4 = greg.ikine(T3,q3); % Numerical solution

% smooth trajectory
% t = [0:0.05:2];
% q = jtraj(q1,q2,t);
% greg.plot(q2);










