function [ qt ] = Greg_IK4( EE,Orient, L1, L2, L3 )
%geometric solution to the IK for 4 DOF

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
    %disp('quadrant 3 or 4');
    theta0 = theta0 + pi;
    theta1 = pi - theta1;
    theta2 = -theta2;
end

theta3 = theta1 - theta2 - Orient;

qt = [theta0, -theta1, theta2,theta3] ;


end

