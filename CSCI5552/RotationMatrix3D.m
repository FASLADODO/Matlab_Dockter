function [Rx, Rz, Ry] = RotationMatrix3D(roll,pitch,yaw)
    Rx = [1,0,0;
        0, cos(roll), -sin(roll);
        0, sin(roll), cos(roll) ];

    Ry = [cos(pitch), 0, sin(pitch);
        0, 1, 0;
        -sin(pitch), 0, cos(pitch) ];

    Rz = [cos(yaw), -sin(yaw), 0;
        sin(yaw), cos(yaw), 0;
        0, 0, 1 ];
end