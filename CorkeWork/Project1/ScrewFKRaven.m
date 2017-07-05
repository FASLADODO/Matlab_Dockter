function [ T ] = ScrewFKRaven( jointvars, jointLengths )
%Compute forward kinematics use matrix exponentials (for Raven)

%joint offset angles
a13 = jointLengths.a13;
a35 = jointLengths.a35;

%make xi vectors
xi_1 = [0;0;0;0;0;1];
xi_3 = [0;0;0;0;sin(a13);cos(a13)];
xi_5 = [0;0;0;0;sin(a13-a35);cos(a13-a35)];

%g_st(0) To get stuff aligned
g_st_zero = [1,0,0,0;
             0,cos(-a13+a35),-sin(-a13+a35),0;
             0,sin(-a13+a35),cos(-a13+a35),0;
             0,0,0,1];


%instead use twistexp() on 6x1 vectors
T = twistexp(xi_1,jointvars(1))*twistexp(xi_3,jointvars(2))*twistexp(xi_5,jointvars(3))*g_st_zero;

end

