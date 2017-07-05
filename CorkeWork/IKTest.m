%% Test IK

EE = [60,40,160];
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

q_out = [theta0, theta1, theta2]*(180/pi)

%%

XD = 60;
YD = 40;
ZD = 160;

h = 72;
a=175;
B = 151;
b=100;
c=100;

A = asin( (sind(B)/b)*a) * (180/pi);
theta3 = 180-A
theta2 = 33.69 + B - 90

