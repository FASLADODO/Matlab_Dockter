function [theta3, theta4, Ax, Ay, Bx, By] = PositionAnalysis4Bar( A1, A2, A3, A4, theta2 )
%Solve position of 4 bar linkage
% Give link lengths Ai and inputs angle beta2 (degrees)
% give beta2 in

%http://facultad.bayamon.inter.edu/elay/mecn4110/ch04%20Position%20Analysis.pdf

Ax = A2 *cosd(theta2);
Ay = A2 *sind(theta2);

%Change of variables

S = (A2^2 - A3^2 + A4^2 - A1^2)/(2*(Ax-A1));
R = (A1 - S)^2 - A4^2;
Q = (2*Ay*(A1-S))/(Ax - A1);
P = (Ay^2)/((Ax-A1)^2) + 1;

By_op = [(-Q + sqrt(Q^2 - 4*P*R)) / (2*P), (-Q - sqrt(Q^2 - 4*P*R)) / (2*P)];
if( isreal(By_op(1)) )
   By = By_op(2);
else
    By = By_op(1);
end
Bx = S - (2*Ay*By)/(2*(Ax-A1));


theta3 = atan2((By - Ay),(Bx - Ax)) * (180/pi);
theta4 = atan2( By , (Bx - A1) ) * (180/pi);

end

