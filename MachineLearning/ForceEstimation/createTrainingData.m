function [input,target,input_true,target_true] = createTrainingData()
%% Building dummy data -- For this data I am creating a output function of force
    % The input and target are the measured values (w/noise)
    % The input_true and target_true are the analytical values (no noise)
% Let's set up a data set
datapts = 60000;                % Number of measurements in the data set
t = 0:1:datapts-1;              % Time vector (dependent on the number of datapts)
freq = 1000;                    % Sampling frequency in Hz
t = t/freq;                     % Convert the time vector to seconds (incremented by 1 ms, i.e.-Sampling at 1 kHz)

% Making a chirp signal
f0 = 0.1;           % Desired starting frequency 0.1 Hz
f1 = 0.5;           % Desired final frequency 5 Hz (10 Hz seemed too high)
T = t(end);         % Duration you want the chirp signal to last
k = (f1-f0)/T;      % Chirpyness (i.e.-how fast your frequency changes)  
phi = 0;            % Phase offset

rate = phi+2*pi*(f0*t+0.5*k*t.^2);        % The rate of your chirp function

noiseLevel = 0.01;                       % This is the percent of noise added to the measurements

%%% IMPORTANT INFO ABOUT CUI ENCODERS BEING SIMULATED !!!
% 2048 counts per revolution w/o quadrature
% 8192 counts per revolution w/ quadrature
COUNTS2RAD = (2*pi)/(2048*4);

% Set up true values for position, velocity, and acceleration
A = 1024;                                               % Amplitude in counts
pos = (A*sin(rate))*COUNTS2RAD;                         % This is the true position in radians
vel = A*2*pi*(f0+k*t).*cos(rate)*COUNTS2RAD;            % This is the true velocity in radians
acc = -A*(2*pi*(f0+k*t)).^2.*sin(rate)*COUNTS2RAD;      % This is the true accel in radians

% Set up measurements for position, velocity, and acceleration
pos_meas = floor(pos/COUNTS2RAD)*COUNTS2RAD;            % Simulating a discretized encoder signal converted to radians
vel_meas = diff(pos_meas);                              % Simple diff to get velocity
vel_meas(end+1) = vel_meas(end);
acc_meas = diff(vel_meas);                              % Simple diff to get accel
acc_meas(end+1) = acc_meas(end);

% Set up torque_in measurements
A1 = 0.025;                                                                        % Amplitude in Newton-meters                            
torque_in = A1*sin(rate);                                                          % This is the torque_in measurements
torque_in_meas = torque_in + noiseLevel*max(torque_in)*randn(size(torque_in)); % Add noise to the true torque

% Set up the torque_out measurements, which is a function of the inputs
% The torque_out is a (made up) function of all the inputs
torque_out = computeTorque(pos,vel,acc,torque_in);
torque_out_meas = torque_out + noiseLevel*max(torque_out)*randn(size(torque_in));

%USE COLUMN MAJOR DOOFUS
% Combine all the inputs into one matrix (X)
input = [pos_meas;vel_meas;acc_meas;torque_in_meas]';
% input = [pos_meas;torque_in_meas]';
target = torque_out_meas';
input_true = [pos;vel;acc;torque_in]';
target_true = torque_out';

end

