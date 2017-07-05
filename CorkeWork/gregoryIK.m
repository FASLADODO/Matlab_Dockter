%Make gregory without wrist
clear all

%offset zero pose
theta_offset = [0,0,pi/2,0,0,0]; %normally would be in DH table

%theta,d,a,alpha (DH convention)
L(1) = Link([0 0 0 -pi/2]); %link 1
L(2) = Link([0 0 100 0]); %link 2
L(3) = Link([0 0 0 pi/2]); %link 3
L(4) = Link([0 100 0 -pi/2]); %link 4
L(5) = Link([0 0 0 pi/2]); %link 5
L(6) = Link([0 0 0 0]); %link 6 

%create robot
greg2 = SerialLink(L, 'name', 'gregory');

% [blank, theta1, theta2, theta3, theta4, theta5,theta6];
qn = [0,0,0,0,0,0] %straight up

greg2.plot(qn + theta_offset)


%% Psuedo IK

% desired pos
EE = [-40,60,160];
Orient = [0,0,0]; %roll,pitch,yaw

%constants mm
L1 = 100;
L2 = 100;


% IK
XE = EE(1);
YE = EE(2);
ZE = EE(3);

%height to wrist
Y_prime = sqrt(XE^2 + YE^2);
h = sqrt( ZE^2 + Y_prime^2);

% base joint 0 rotation
theta0 = atan2(YE,XE);

if( XE^2 + YE^2 + ZE^2 >= (L1 + L2)^2 )
    %If outside of reachable workspace, map back
    disp('out of workspace');
    theta1 = atan2(abs(ZE),abs(Y_prime));
    theta2 = 0;
    
else
    % joint 1 and joint 2
    phi = acos( (L1^2 + L2^2 - h^2) / (2*L1*L2) );
    theta2 = pi - phi;
    alpha = atan2(abs(ZE),Y_prime);
    beta = asin( (L2/h) * sin(phi) );
    theta1 = alpha + beta;
end
    
if (theta0 < 0)
    disp('quadrant 3 or 4');
    theta0 = theta0 + pi;
    theta1 = pi - theta1;
    theta2 = -theta2;
end

theta3 = 0; %ignore yaw
theta4 = theta1 - theta2 - Orient(2);
theta5 = Orient(1);

qt = [theta0, -theta1, theta2,theta3,theta4,theta5] 


greg2.plot(qt + theta_offset)

T0 = greg2.fkine(qt+theta_offset)
