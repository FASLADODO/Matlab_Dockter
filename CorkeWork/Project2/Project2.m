%%Project 1
%Rod Dockter
%Using the Gregory Arm (simplified to 4 DOF) 
%Base Swivel, Shoulder, Elbow, and wrist pitch


%We will optimize for a welding job, assuming the welding gun is affixed at
%the end of the wrist.

%Assume that the optimal welding speed is 0.005 m/s (via the internet)

%We also assume that the welding surface is always a plane so that position
%of the welding gun should be straight down onto plane

%We will optimaize to be able to weld in a 5cm x 5cm region

%QUESTION A
%I will assume that we are using the Tower Pro MG995 servos found in Gregory
%these servos have an angular velocity of:
%0.2 sec/60 deg = 300 deg/sec = 5.236 rad/s
% Therefore to achieve 1 m/s the effective radius would need to be v = r*w
% in other words radius = v / w = (1 m/s)/(5.236 rad/s) = 0.191 m

%----------------------------------------------------------------------------------%
%% Just for plotting

clear all

L1 = 100;
L2 = 100;
L3 = 20;

%offset zero pose
theta_offset = [0,-pi/2,0,pi/2]; %normally would be in DH table

%theta,d,a,alpha (DH convention)
L(1) = Link([0 0 0 -pi/2]); %link 1
L(2) = Link([0 0 L1 0]); %link 2
L(3) = Link([0 0 L2 0]); %link 3
L(4) = Link([0 0 0 pi/2]); %link 4

%create robot
greg4 = SerialLink(L, 'name', 'gregory4');

greg4.tool = transl(0,0,L3);

% [blank, theta1, theta2, theta3, theta4, theta5,theta6];
%qn = [0,-pi/4,3*pi/4,0]; %straight up
qn = [0,0,0,0]; %straight up

greg4.plot(qn )

title('Welding Position')

%----------------------------------------------------------------------------------%
%% Getting jacobian

clear all

%We have 3 link lengths we can vary (Upper arm, lower arm, and wrist
%offset)
syms L1 L2 L3 q1 q2 q3 q4

%offset zero pose
theta_offset = [0,-pi/2,0,pi/2]; %normally would be in DH table

%theta,d,a,alpha (DH convention)
L(1) = Link([0 0 0 -pi/2]); %link 1
L(2) = Link([0 0 L1 0]); %link 2
L(3) = Link([0 0 L2 0]); %link 3
L(4) = Link([0 0 0 pi/2]); %link 4

%create robot
greg4 = SerialLink(L, 'name', 'gregory4')

%Wrist offset
greg4.tool = transl(0,0,L3);

% [blank, theta1, theta2, theta3, theta4, theta5,theta6];
qt = [q1,q2,q3,q4]; %straight up

%get symbolic FK and jacobian
FKsym   = greg4.fkine ( qt )                    
Jsym    = greg4.jacob0( qt) 
Jsym    = Jsym([1:3,5], 1:4)       % we only want x,y,z,pitch     
vpa(det(Jsym),2)                        

%----------------------------------------------------------------------------------%
%%
% Singular Values wont return an answer for 4 DOF symbolic :-(
% disp('Get singular values (symbolic!):')
% sigma = simplify( eig(Jsym.' * Jsym));
% s1 = sigma(1) ;  % get the lowest singular value
% pretty(sigma)

%% Now find optimal manipulability region for generic lengths
% Initial results indicated that manipulability was greatest when base
% swivel was zero q1 = 0

%Additionaly the wrist manipulability change was minimal for most angle
%Therefore we will focus on the shoudler and elbow joints

L1 = 100;
L2 = 100;
L3 = 20;

%theta,d,a,alpha (DH convention)
L(1) = Link([0 0 0 -pi/2]); %link 1
L(2) = Link([0 0 L1 0]); %link 2
L(3) = Link([0 0 L2 0]); %link 3
L(4) = Link([0 0 0 pi/2]); %link 4

%create robot
greg4 = SerialLink(L, 'name', 'gregory4')

%Wrist offset
greg4.tool = transl(0,0,L3);

Jdet = [];
Manip = [];
inder = 1;
for ii = [0:.4:2*pi]
    for jj = [0:.4:2*pi]
        for kk = [0:.4:2*pi]
            for ll = 4.8
                J = greg4.jacob0([ii,jj,kk,ll]);
                J = J([1:3,5], 1:4); %just xyz,pitch
                S = svd(J);
                T = greg4.fkine([ii,jj,kk,ll]);
                x = T(1,end);
                y = T(2,end);
                z = T(3,end);
                Manip(inder,:) = [ii,jj,kk,ll,x,y,z,min(S),max(S)];
                %Manip(ii,jj,kk,ll)=min(S);
                inder = inder + 1;
            end
        end
    end
end
max(Manip(:,8))


%----------------------------------------------------------------------------------%
%% Plot the manipulability measure min(sigma) for all q1 and q2 angles
figure(1)

subplot(2,1,1)
plot3(Manip(:,2),Manip(:,3), Manip(:,8),        'rx'); hold on
%plot3(Manip(:,2),Manip(:,3), Manip(:,9),        'gx'); hold on %%max
plot3(Manip(:,2),Manip(:,3), Manip(:,8)./Manip(:,9),'b.'); hold off

subplot(2,1,2)
plot3(Manip(:,5),Manip(:,6), Manip(:,8),       'rx'); hold on
%plot3(Manip(:,5),Manip(:,6), Manip(:,9),        'gx'); hold on %%max
plot3(Manip(:,5),Manip(:,6), Manip(:,8)./Manip(:,9),'b.'); hold off

subplot(2,1,1)
grid on
xlabel('q1')
ylabel('q2')
zlabel('Manipulability')
legend('min','min/max')

subplot(2,1,2)
grid on
xlabel('x')
ylabel('y')
zlabel('Manipulability')
legend('min','min/max')

%Results show that the region near the origin has the most manipulability
% 0<x<100 -25<y<25 z = 0
%ie  70 < q2 < 90 and 150 < q3 < 180
%Additionally the q1 = 0 direction has the most manipulability of all
%directions. This makes sense since the wrist is perpendicular to that
%direction


%----------------------------------------------------------------------------------%
%% Find optimal L1, L2, L3 using min(sigma)
% Assume for gregory that the optimal manipulability happens at
%q_op = [0,-pi/4,3*pi/4,0]

clear all

syms L1 L2 L3 q1 q2 q3 q4

%offset zero pose
theta_offset = [0,-pi/2,0,pi/2]; %normally would be in DH table

%theta,d,a,alpha (DH convention)
L(1) = Link([0 0 0 -pi/2]); %link 1
L(2) = Link([0 0 L1 0]); %link 2
L(3) = Link([0 0 L2 0]); %link 3
L(4) = Link([0 0 0 pi/2]); %link 4

%create robot
greg4 = SerialLink(L, 'name', 'gregory4')

%Wrist offset
greg4.tool = transl(0,0,L3);

% [blank, theta1, theta2, theta3, theta4, theta5,theta6];
qt = [q1,q2,q3,q4]; %straight up %straight up

%get symbolic FK and jacobian                 
Jsym    = greg4.jacob0( qt) 
Jsym    = Jsym([1:3,5], 1:4)       % we only want x,y,z,pitch     
det(Jsym)  

%sub in the values
q1 = 0;
q2 = 3*pi/4;
q3 = -pi/4;
q4 = pi/2;
JS = eval(Jsym);   % Sub in qi into J
JS = vpa( JS,2)  % approximate to decimal

% Find manipulability dependence on L1 and L2 lengths


%loop through and recompute sigma at each L1 L2 L3
manip_length = [];
nnd = 1;
for L1 = 0:5:100
    for L2 = 0:5:100
        for L3 = 0:5:100
            JN = eval(JS);   % substitute in L1 L2 L3
            sig = svd(JN);
            manip_length(nnd,:) = [L1, L2, L3, min(sig)];
            nnd = nnd +1;
        end
    end
end

[max_sig, index] = max(manip_length(:,4));

fprintf('Highest manipulability occurs at L1,L2,L3 = \n');
manip_length(index,1:3)
fprintf('manblty= \n')
max_sig


%Note: As expected, the highest manipulability comes from L1=L2=100 L3=0 (the max
% length)
% The manipulability measure here is min 0.99



%----------------------------------------------------------------------------------%
%% Plot all the link lengths vs manipulability

figure(21)
plot3(manip_length(:,1),manip_length(:,2),manip_length(:,4),'b+');
xlabel('L1 (mm)')
ylabel('L2 (mm)')
zlabel('sig val')
figure(22)
plot3(manip_length(:,1),manip_length(:,3),manip_length(:,4),'g+');
xlabel('L1 (mm)')
ylabel('L3 (mm)')
zlabel('sig val')
figure(23)
plot3(manip_length(:,2),manip_length(:,3),manip_length(:,4),'r+');
xlabel('L2 (mm)')
ylabel('L3 (mm)')
zlabel('sig val')



%----------------------------------------------------------------------------------%
%% QUESTION B 
% We will use the values for L1, L2 and L3 found above that maximize the
% manipulability

% Next we want to perform welding in a 5cmx5cm (50mm x 50mm) square.
% In order to perform this action in the region of optimal manipulability
% we will perform this action in the region
% 0<x<50 -25<y<25 z = 0 (mm)
% Additionally the weld needs to point straight down ie
% pitch = -90

%In order to prove that the required speed ( 0.005 m/s) is achievable
%throughout this region with the chosen link lengths, the inverse
%kinematics will be solved for a grid of points in the desired region.
% With the angles found from IK the jacobian will be solved for.
% Then with that jacobian, the cartesian velocity will be computed using
% the max joint velocities in order to ensure that the the max velocity is
% always greater than the required velocity

clear all

%link lengths 
L1 = 100; %(100 mm)
L2 = 100; %(100 mm)
L3 = 10; %note L3=0 because that gave the highest manipulability

% joint velocity max vector (5.236 rad/s)
q_dot = [5.236;5.236;5.236;5.236];

%theta,d,a,alpha (DH convention)
L(1) = Link([0 0 0 -pi/2]); %link 1
L(2) = Link([0 0 L1 0]); %link 2
L(3) = Link([0 0 L2 0]); %link 3
L(4) = Link([0 0 0 pi/2]); %link 4

%create robot
greg4 = SerialLink(L, 'name', 'gregory4')

%Wrist offset
greg4.tool = transl(0,0,L3);

%iterate through all the xyz points (in m
step = 2;
X = 50:step:100;
Y = -25:step:25;
Z = 0;
pitch = pi/2; %straight down

%store all positions and velocities
Velocity_Mat = [];
indv = 1;

% loop through discrete positions in workspace
for xx = X
    for yy = Y
        q_out = Greg_IK4( [xx,yy,Z] ,-pi,L1,L2,L3);
        %get jacobian and FK
        FK   = greg4.fkine ( q_out );
        J_mat    = greg4.jacob0( q_out ); %(in mm units)
        %get cartesian velocity and pack into mat
        %cart_vel = J_mat*q_dot./1000; %to get into m/s
        svder = svd(J_mat([1:3,5],:));
        Velocity_Mat(indv,:) = [FK(1,end),FK(2,end),FK(3,end),q_out,min(svder)/1000];
        indv = indv + 1;
    end
end

%Now plot all the points in the grid along with the minimum value of the
%cartesian velocity
plot3(Velocity_Mat(:,1),Velocity_Mat(:,2),Velocity_Mat(:,8),'bx')
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('min velocity (m / s)')
title('Minimum cartesian velocity in workspace (0.005 m/s required)')


%This results in an overall minimum velocity of 0.05 m/s in the workspace
%This exceeds our minimum requirement



%----------------------------------------------------------------------------------%
%% QUESTION C 
% The graph clearly indicates that the minimum velocity is always 
% above the required value of 0.005 m/s
% This means that we can choose a smaller set of values for L1 and L2

% Based on the manip_length() matrix it appears that L1=L2=60 still results
% in a manipulability of 0.996 which is close to the maximum.
% By inspection, the minimum link length to use this workspace is
% sqrt(100^2+25^2)/2 = 52 mm, therefore L1=L2=60 is still value

clear all

%New link lengths 
L1 = 60;
L2 = 60;
L3 = 0; %note L3=0 because that gave the highest manipulability

% joint velocity max vector (5.236 rad/s)
q_dot = [5.236;5.236;5.236;5.236];

%theta,d,a,alpha (DH convention)
L(1) = Link([0 0 0 -pi/2]); %link 1
L(2) = Link([0 0 L1 0]); %link 2
L(3) = Link([0 0 L2 0]); %link 3
L(4) = Link([0 0 0 pi/2]); %link 4

%create robot
greg4 = SerialLink(L, 'name', 'gregory4')

%Wrist offset
greg4.tool = transl(0,0,L3);

%iterate through all the xyz points (in m
step = 2;
X = 50:step:100;
Y = -25:step:25;
Z = 0;
pitch = pi/2; %straight down

%store all positions and velocities
Velocity_Mat = [];
indv = 1;

% loop through discrete positions in workspace
for xx = X
    for yy = Y
        q_out = Greg_IK4( [xx,yy,Z] ,-pi,L1,L2,L3);
        %get jacobian and FK
        FK   = greg4.fkine ( q_out );
        J_mat    = greg4.jacob0( q_out ); %(in mm units)
        %get cartesian velocity and pack into mat
        svder = svd(J_mat(1:3,1:3));
        Velocity_Mat(indv,:) = [FK(1,end),FK(2,end),FK(3,end),q_out,min(svder)/1000];
        indv = indv + 1;
    end
end

%Now plot all the points in the grid along with the minimum value of the
%cartesian velocity
plot3(Velocity_Mat(:,1),Velocity_Mat(:,2),Velocity_Mat(:,8),'bx')
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('min velocity (m / s)')
title('Minimum cartesian velocity in workspace (with 60 mm link lengths)')

% All locations in the grid still result in a sufficient cartesian velocity
% Minimum velocity in workspace is now 0.027 m/s
% Therefore we find that the L1 = L2 = 60mm, L3 = 0 link combination
% approximates the smallest link lengths that still provide sufficient
% speed for the welding application

% It should be noted that L3 = 0 is slightly unrealistic since this implies
% the welding device is concentric with the wrist but oh well.


%There is likely not a closed form solution for the minimum link lengths
%The symbolic SVD for the jacobian will be a nonlinear equation with 7
%variables for which there is rarely a closed form solution.