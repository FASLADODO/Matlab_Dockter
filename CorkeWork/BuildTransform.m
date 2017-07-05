function [ T ] = BuildTransform( tx,ty,tz,roll,pitch,yaw )
%Build 4x4 transformation matrix from position and roll, pitch, yaw
% angles in radians

% Use
% T = BuildTransform( 10,5,7,0.1,0.2,0.3 );

T = eye(4);
R = eul2r(roll, pitch, yaw);
P = [tx;ty;tz];
T(1:3,1:3) = R;
T(1:3,4) = P;

end

