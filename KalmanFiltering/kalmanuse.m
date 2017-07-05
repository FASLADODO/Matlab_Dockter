% 1D motion with no input
clear s
s.x = 0;
s.A = 1;
% Define a process noise (stdev) of 2 volts as the car operates:
s.Q = 0.8^2; % variance, hence stdev^2
% Define the voltimeter to measure the voltage itself:
s.H = 1;
% Define a measurement error (stdev) of 2 volts:
s.R = 0.8^2; % variance, hence stdev^2
% Do not define any system input (control) functions:
s.B = 0;
s.u = 0;
% Do not specify an initial state:
s.x = 0;
s.P = nan;
% Generate random voltages and watch the filter operate.
tru=[]; % truth voltage
for t=1:100
   tru(end+1) = randn*2+12;
   s(end).z = tru(end) + randn*2; % create a measurement
   s(end+1)=kalmanf(s(end)); % perform a Kalman filter iteration
end
figure
hold on
grid on
for k = 1:100
	% plot measurement data:
    hz=plot([s(k).z],0,'r.');
    % plot a-posteriori state estimates:
    hk=plot([s(k).x],0,'b.');
    ht=plot(tru(k),0,'g.');
	axis equal
    pause(1);
	M(k) = getframe;
    clf;
end


movie(M,30,1)
