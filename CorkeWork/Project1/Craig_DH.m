function [ T ] = Craig_DH( theta, d, a, alpha )
%Implements a transformation matrix, using craig convention
%This convention is dumb, Spong 4 lyfe

T = [          cos(theta),           -sin(theta),           0,             a;
    sin(theta)*cos(alpha), cos(theta)*cos(alpha), -sin(alpha), -d*sin(alpha);
    sin(theta)*sin(alpha), cos(theta)*sin(alpha),  cos(alpha),  d*cos(alpha);
                        0,                     0,           0,            1];

end

