function [ T ] = ScrewFKRhino( jointvars )
%Compute forward kinematics use matrix exponentials (For Rhino)

%joint lengths
d1 = 195;
a2 = 170; 
a3 = 170;
a4 = 1; 
d5 = 125;

%make xi vectors
xi_1 = [0;0;0;0;0;1];
xi_2 = [-d1;0;0;0;1;0];
xi_3 = [-(a2+d1);0;0;0;1;0];
xi_4 = [-(d1+a2+a3);0;0;0;1;0];
xi_5 = [0;a4;0;0;0;1];

%g_st(0)
g_st_zero = [1,0,0,-a4;
            0,1,0,0;
            0,0,1,d1+a2+a3+d5;
            0,0,0,1];

%instead use twistexp() on 6x1 vectors
T = twistexp(xi_1,jointvars(1))*twistexp(xi_2,jointvars(2))*twistexp(xi_3,jointvars(3))*twistexp(xi_4,jointvars(4))*twistexp(xi_5,jointvars(5))*g_st_zero;


end

