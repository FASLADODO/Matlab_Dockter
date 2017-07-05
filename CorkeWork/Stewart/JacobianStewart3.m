function [A] = JacobianStewart3(r_p,r_base,A1,A2)

syms x y z ;

% Setup length
L  = r_base ; 				% Base
l  = r_p	;			% Platform
l1 = A1 ;				% Link1
l2 = A2 ;				% Link2

% Motor1
	a1 = x^2 + y^2 + z^2 + 2*y*l/sqrt(3) - 2*y*L/sqrt(3) + l^2/3 - 2*l*L/3 + L^2/3 + l1^2 - l2^2 ;
	b1 = 2*y*l1 + 2*l*l1/sqrt(3) - 2*l1*L/sqrt(3) ;
	c1 = 2*z*l1 ;
	% Select Configuration
    %q1 = atan2(c1,b1) + acos(a1/sqrt(b1^2+c1^2)) ;		% Upper Configuration
    q1 = atan(c1/b1) - acos(a1/sqrt(b1^2+c1^2)) ;		% Lower Configuration
   
% Motor2
	a2 = x^2 + y^2 + z^2 + x*l - x*L - y*l/sqrt(3) + y*L/sqrt(3) + l^2/3 - 2*l*L/3 + L^2/3 + l1^2 - l2^2 ;
	b2 = sqrt(3)*x*l1 + 2*l*l1/sqrt(3) - 2*L*l1/sqrt(3) - y*l1 ;
    c2 = 2*z*l1 ;
    % Select Configuration
    %q2 = atan2(c2,b2) + acos(a2/sqrt(b2^2+c2^2)) ;		% Upper Configuration
    q2 = atan(c2/b2) - acos(a2/sqrt(b2^2+c2^2)) ;		% Lower Configuration
   
% Motor3
    a3 = x^2 + y^2 +z^2 - x*l + x*L - y*l/sqrt(3) + y*L/sqrt(3) + l^2/3 - 2*l*L/3 + L^2/3 + l1^2 - l2^2 ;
    b3 = - sqrt(3)*x*l1 - y*l1 + 2*l*l1/sqrt(3) - 2*L*l1/sqrt(3) ;
    c3 = 2*z*l1 ;
    % Select Configuration
    %q3 = atan2(c3,b3) + acos(a3/sqrt(b3^2+c3^2)) ;		% Upper Configuration
    q3 = atan(c3/b3) - acos(a3/sqrt(b3^2+c3^2)) ;		% Lower Configuration

A = jacobian([q1 ; q2 ; q3] , [x y z]) ;