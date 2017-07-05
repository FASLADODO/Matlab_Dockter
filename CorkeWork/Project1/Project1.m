%%%%% Project 1   %%%%%
%%%%% ME 8287     %%%%%
%%%%% Rod Dockter %%%%%
%%%%% Code        %%%%%


%% Cross products for xi vectors

syms d1 a2 a3 a4 d5 s13 c13

q = [a4;0;d1+a2+a3];
w = [0;0;1];

csi = cross(w,q)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RHINO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Question A, part c)

% Implementing Rhino in Matrix exponentials

clear all

%this is one indexed (ie 1-5) not (0-4)
syms d1 a2 a3 a4 d5 real
syms q_1 q_2 q_3 q_4 q_5 real
syms s_1 s_2 s_3 s_4 s_5 real
syms c_1 c_2 c_3 c_4 c_5 real
% theta1=0;theta2=0;theta3=0;theta4=0;theta5=0;


%make vectors
xi_1 = [0;0;0;0;0;1];
xi_2 = [-d1;0;0;0;1;0];
xi_3 = [-(a2+d1);0;0;0;1;0];
xi_4 = [-(d1+a2+a3);0;0;0;1;0];
xi_5 = [0;a4;0;0;0;1];

%create skew matrices (NOT USED)
xi_hat_1 = twist(xi_1);
xi_hat_2 = twist(xi_2);
xi_hat_3 = twist(xi_3);
xi_hat_4 = twist(xi_4);
xi_hat_5 = twist(xi_5);

%g_st(0)
g_st_zero = [1,0,0,-a4;
            0,1,0,0;
            0,0,1,d1+a2+a3+d5;
            0,0,0,1];
        
% using twistexp() on twist() matrices doesnt work
% g_st = twistexp(xi_hat_1,theta1)*twistexp(xi_hat_2,theta2)*twistexp(xi_hat_3,theta3)*twistexp(xi_hat_4,theta4)*twistexp(xi_hat_5,theta5)*g_st_zero

%instead use twistexp() on 6x1 vectors
g_st = twistexp(xi_1,q_1)*twistexp(xi_2,q_2)*twistexp(xi_3,q_3)*twistexp(xi_4,q_4)*twistexp(xi_5,q_5)*g_st_zero;

%simplifying everything for better print out cos(q_1) = c_1
g_st = simplify(g_st);
g_st = subs(g_st,sin(q_1),s_1);
g_st = subs(g_st,sin(q_2),s_2);
g_st = subs(g_st,sin(q_3),s_3);
g_st = subs(g_st,sin(q_4),s_4);
g_st = subs(g_st,sin(q_5),s_5);
g_st = subs(g_st,cos(q_1),c_1);
g_st = subs(g_st,cos(q_2),c_2);
g_st = subs(g_st,cos(q_3),c_3);
g_st = subs(g_st,cos(q_4),c_4);
g_st = subs(g_st,cos(q_5),c_5);

fprintf('Rhino Robot Transform \n \n T = \n \n');
pretty(g_st) % See print out

%% Implementing Rhino in Corke (just to plot and check)

clear all

%joint lengths
d1 = 195;
a2 = 170; 
a3 = 170;
a4 = 1; 
d5 = 125;

%offset zero pose, so it makes robot straight up in air, matches screw
%formulation
theta_offset = [0,-pi/2,0,-pi/2,pi]; %normally would be in DH table

%theta,d,a,alpha (DH convention)
L(1) = Link([0 d1 0 -pi/2]); %link 1
L(2) = Link([0 0 a2 0]); %link 2
L(3) = Link([0 0 a3 0]); %link 3
L(4) = Link([0 0 a4 -pi/2]); %link 4
L(5) = Link([0 d5 0 0]); %link 5

%create robot
Rhino = SerialLink(L, 'name', 'Rhino');

% [blank, theta1, theta2, theta3, theta4, theta5,theta6];
q_in = [0,0,0,0,0]; %straight up

%plot robot for pretty sake
Rhino.plot(q_in+theta_offset)

%get transformation matrix for given joint variables
Rhino.fkine(q_in+theta_offset)

%% Question A, part d)
% Should result in small errors (~2^-15)

%Comparing Corke and Matrix Exp for Rhino

max = 2*pi; %0-360 degrees
nn = 1000;
%create random joint variables, combine into matrix
q1 = max.*rand(nn,1);
q2 = max.*rand(nn,1);
q3 = max.*rand(nn,1);
q4 = max.*rand(nn,1);
q5 = max.*rand(nn,1);

q_t = [q1,q2,q3,q4,q5];

%loop through and compare random joint configuration for DH vs exponential
%formulation
for ii = 1:nn
    %using custom functions for FK calculations
    T_exp = ScrewFKRhino(q_t(ii,:));
    T_DH = DHFKRhino(q_t(ii,:));
    
    %compute sum of differences in transformation matrices
    diff(ii) = sum(sum(abs(T_exp - T_DH)));
end

%plot differences
plot(1:nn,diff);
xlabel('Iteration #')
ylabel('Error (mixed units)')
title('Error between DH and Exponential Formulation (Rhino)')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Raven %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Question A, part c)

% Implementing Raven in Matrix exponentials

clear all

%this is one indexed (ie 1,3,5) not (0-2) because of Lum notation
syms a13 a35 real
syms c_13 c_35 c_135 s_13 s_35 s_135 real
syms q_1 q_3 q_5  real
syms s_1 s_3 s_5 real
syms c_1 c_3 c_5 real

%q = [0,0,0];
q = [q_1, q_3, q_5];

%make vectors
xi_1 = [0;0;0;0;0;1];
xi_3 = [0;0;0;0;sin(a13);cos(a13)];
xi_5 = [0;0;0;0;sin(a13-a35);cos(a13-a35)];

%g_st(0)
g_st_zero = [1,0,0,0;
             0,cos(-a13+a35),-sin(-a13+a35),0;
             0,sin(-a13+a35),cos(-a13+a35),0;
             0,0,0,1];
        
%instead use twistexp() on 6x1 vectors
g_st = twistexp(xi_1,q(1))*twistexp(xi_3,q(2))*twistexp(xi_5,q(3))*g_st_zero;

%simplifying everything for better print out cos(q_1) = c_1
g_st = simplify(g_st);
g_st = subs(g_st,sin(q_1),s_1);
g_st = subs(g_st,sin(q_3),s_3);
g_st = subs(g_st,sin(q_5),s_5);
g_st = subs(g_st,cos(q_1),c_1);
g_st = subs(g_st,cos(q_3),c_3);
g_st = subs(g_st,cos(q_5),c_5);

g_st = subs(g_st,sin(a13),s_13);
g_st = subs(g_st,sin(a35),s_35);
g_st = subs(g_st,cos(a13),c_13);
g_st = subs(g_st,cos(a35),c_35);

g_st = subs(g_st,sin(a13 + a35),s_135);
g_st = subs(g_st,cos(a13 + a35),c_135);

fprintf('Raven Robot Transform \n \n T = \n \n');
pretty(g_st) % See print out

%% Implementing Raven in Corke, DOES NOT WORK!, LUM USES CRAIG!
% What a n00b!

clear all

%joint offset angles
a13 = 45*(pi/180);
a35 = 90*(pi/180); 

%theta,d,a,alpha (DH convention)
L(1) = Link([0 0 0 a13]); %link 1
L(2) = Link([0 0 0 a35]); %link 2
L(3) = Link([0 0 0 0]); %link 3

%create robot
Raven = SerialLink(L, 'name', 'Raven');

% [blank, theta1, theta2, theta3, theta4, theta5,theta6];
q_in = [0,0,0]; %straight up

%plot robot for pretty sake
%Raven.plot(q_in)

%get transformation matrix for given joint variables
Raven.fkine(q_in)


%% Implementing Raven in hacky Craig DH method

%Test values for parameters
jointLengths.a13 = 45*(pi/180);
jointLengths.a35 = 90*(pi/180); 

% Simple pose
q_in = [0,0,0]; %straight up

%Get transformation matrix for FK
T_fk = DHFKRaven(q_in, jointLengths)



%% Question A, part d)
% Should result in small errors (~2^-15)

%Comparing Corke and Matrix Exp FK for Raven
jointLengths.a13 = pi/5;
jointLengths.a35 = pi/3; 

max = 2*pi; %0-360 degrees
nn = 1000;
%create random joint variables, combine into matrix
q1 = max.*rand(nn,1);
q3 = max.*rand(nn,1);
q5 = max.*rand(nn,1);
q_t = [q1,q3,q5];

%loop through and compare random joint configuration for DH vs exponential
%formulation
for ii = 1:nn
    %using custom functions for FK calculations
    T_exp = ScrewFKRaven(q_t(ii,:),jointLengths);
    T_DH = DHFKRaven(q_t(ii,:),jointLengths);
    
    %compute sum of differences in transformation matrices
    diff(ii) = sum(sum(abs(T_exp - T_DH)));
end

%plot differences
plot(1:nn,diff);
xlabel('Iteration #')
ylabel('Error (mixed units)')
title('Error between DH and Exponential Formulation (Raven)')






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Inverse Kinematics %%%%%%%%%%%%%%%%%%%%

%% Question B, part a)

%Solving the Gregory Arm using the corke IK method, this is for comparison with
%the geometric method by hand

clear all

%joint parameters (symbolic)
syms L1 L2 real

%offset zero pose
theta_offset = [0,-pi/2,pi/2,0,0,0]; %normally would be in DH table

%theta,d,a,alpha (DH convention)
L(1) = Link([0 0 0 -pi/2]); %link 1
L(2) = Link([0 0 L1 0]); %link 2
L(3) = Link([0 0 0 pi/2]); %link 3
L(4) = Link([0 L2 0 -pi/2]); %link 4
L(5) = Link([0 0 0 pi/2]); %link 5
L(6) = Link([0 0 0 0]); %link 6 

%create robot
Greg = SerialLink(L, 'name', 'Greg');

Greg.tool = transl(0,0,0); %end effector offset

% [blank, theta1, theta2, theta3, theta4, theta5,theta6];
qn = [0,0,0,0,0,0]; %straight up
qw = [0,0,pi/2,0,-pi/2,0]; %elbow out, wrist up

%Plot particular positions
%Greg.plot(qn + theta_offset)

% Symbolic IK solution
q_out = Greg.ikine_sym(6)

fprintf('Question B, part b), Symbolic Gregory IK Solution: \n \n');
% solution is legit, solution 1 is not? Who knows?
s1 = q_out{1}(2)  % is one solution
s2 = q_out{2}(2)
s3 = q_out{3}(2)
s4 = q_out{4}(2)  % is one solution
s5 = q_out{5}(2)
s6 = q_out{6}(2)

%See PDF print out

%% Testing inverse kinematics solution
% Implementing gregory in corke and then using a pose, get the transform, 
% plug it back in for inverse kinematics compare joint angles

clear all

L1 = 100;
L2 = 100;
jointLengths.L1 = L1;
jointLengths.L2 = L2;

%offset zero pose
theta_offset = [0,-pi/2,pi/2,0,0,0]; %normally would be in DH table

%theta,d,a,alpha (DH convention)
L(1) = Link([0 0 0 -pi/2]); %link 1
L(2) = Link([0 0 L1 0]); %link 2
L(3) = Link([0 0 0 pi/2]); %link 3
L(4) = Link([0 L2 0 -pi/2]); %link 4
L(5) = Link([0 0 0 pi/2]); %link 5
L(6) = Link([0 0 0 0]); %link 6 

%create robot
Greg = SerialLink(L, 'name', 'Greg');

Greg.tool = transl(0,0,0); %end effector offset

% test configuration
qw = [0,0,pi/3,0,-pi/2,0]; %elbow out, wrist up

%Get FK transformation for particular joint angles
fprintf('Forward Kinematics T matrix: \n')
T_out = Greg.fkine(qw+theta_offset)

%Using symbolic IK solution, get back joint angles for that transformation
q_out = GregoryIK_sym( T_out, jointLengths );

fprintf('\n compared with the Inverse -> forward kinematics T: \n')

%Get forward T matrix from inverse solution joing angles
T_test = Greg.fkine(q_out)

fprintf(' \n One is just the elbow up vs down version \n')

%Greg.plot(qw + theta_offset);
Greg.plot(q_out);
title('Inverse Kinematics Pose')
