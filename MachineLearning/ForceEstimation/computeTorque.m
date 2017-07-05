function [ torque_out ] = computeTorque(pos,vel,acc,torque_in)
% This is a made up function to compute the torque at the end effector

% Set up parameters
m = 0.005;          % This is the mass (kg)
b1 = 0.001;         % This is the damping coefficient (N-m*s)
k = 0.2;            % This is the spring coefficient (N-m)

term1 = m*acc.^2;
term2 = b1*vel+sign(vel).*vel;
term3 = k*pos+k*pos.^2;
term4 = 0.5*torque_in;
torque_out = m*acc.^2+b1*vel+sign(vel).*vel+k*pos+k*pos.^2-0.5*torque_in;


Apower = [abs(term1) ./ abs(torque_out); abs(term2) ./ abs(torque_out); abs(term3) ./ abs(torque_out); abs(term4) ./ abs(torque_out)]';

figure
area(Apower)
legend('acc','vel','pos','T')
ylim([0 2])
end

