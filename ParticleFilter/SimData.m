%make data

clear all
close all
clc


%Taken from Astrini thesis, typical alpha, beta and d vals
tissueOptions = 4;
alphaSi = [492.91,489.36,0.0046,356.6571];
betaSi = [14.55,13.79,20.59,11.41];
dSi = [3.1,3.3,1.3,0.5];
typeSi = {'liver','sm-bowel','bladder','gallbladder'};

actualIndex = 1;

% Dynamic Parameters
mYu = .005; %0.005-Yu;
kYu = 0; %0.005-Yu;

dYu = dSi(actualIndex);
alphaYu = alphaSi(actualIndex);
betaYu = betaSi(actualIndex);

%store runs
store_x1 = [];
store_xdot1 = []
store_xdotdot1 = []
store_x2 = [];
store_xdot2 = []
store_xdotdot2 = []

noise = 0.05;

runs = 10;
for ii = 1:runs

    Fgrasp = 2 ;% Hz;
    A = 2; % 5 / (56.3 ) *1000^2; % grasp amplitude (...stress = Applied Newtons/JawArea)A = 2 ; % [newtons]
    % Input to the system, u = force
    tend = 0.25;
    T = 0.001; % sampling period is fronm 1KHz
    t = 0:0.001:tend;
    %    i = int8(1+(1000*t));
    input.time = t;
    input.signals.values = A*  sin(  2*pi*Fgrasp*t );
        %input.signals.values(i) = 50;

    input.time = [input.time]';
    input.signals.values = [input.signals.values]';
    input.signals.dimensions = 1;

    % SIMULINK
    fprintf('Running simulation ...')
    %sim('Astrini_YuNonLinearTissueWithStates.mdl');
    sim('Astrini_YuNonLinearTissueWithStatesWithNoise.mdl');
    fprintf(' DONE\n')

    % no noise
    u = force(:,1); %force.Data ;
    x = state(:,3); %state.Data(:,3) ; 
    xdotsim = state(:,2); %state.Data(:,3) ; 
    xdotdotsim = state(:,1); %state.Data(:,3) ; 
    state = state;

    u = u + (rand(size(u))-.5)*noise*A;    % add noise 
    x = x + (rand(size(x))-.5)*noise*max(x); 
    xdotsim = xdotsim + (rand(size(xdotsim))-.5)*noise*max(xdotsim); 
    xdotdotsim = xdotdotsim + (rand(size(xdotdotsim))-.5)*noise*max(xdotdotsim); 
    store_x1 = [store_x1;x(length(x)/2:end)];
    store_xdot1 = [store_xdot1; xdotsim(length(x)/2:end)];
    store_xdotdot1 = [store_xdotdot1; xdotdotsim(length(x)/2:end)];
end

%change parameters
actualIndex = 2;

dYu = dSi(actualIndex);
alphaYu = alphaSi(actualIndex);
betaYu = betaSi(actualIndex);

for ii = 1:runs

    Fgrasp = 2 ;% Hz;
    A = 2; % 5 / (56.3 ) *1000^2; % grasp amplitude (...stress = Applied Newtons/JawArea)A = 2 ; % [newtons]
    % Input to the system, u = force
    tend = 0.25;
    T = 0.001; % sampling period is fronm 1KHz
    t = 0:0.001:tend;
    %    i = int8(1+(1000*t));
    input.time = t;
    input.signals.values = A*  sin(  2*pi*Fgrasp*t );
        %input.signals.values(i) = 50;

    input.time = [input.time]';
    input.signals.values = [input.signals.values]';
    input.signals.dimensions = 1;

    % SIMULINK
    fprintf('Running simulation ...')
    %sim('Astrini_YuNonLinearTissueWithStates.mdl');
    sim('Astrini_YuNonLinearTissueWithStatesWithNoise.mdl');
    fprintf(' DONE\n')

    % no noise
    u = force(:,1); %force.Data ;
    x = state(:,3); %state.Data(:,3) ; 
    xdotsim = state(:,2); %state.Data(:,3) ; 
    xdotdotsim = state(:,1); %state.Data(:,3) ; 
    state = state;

    u = u + (rand(size(u))-.5)*noise*A;    % add noise 
    x = x + (rand(size(x))-.5)*noise*max(x); 
    xdotsim = xdotsim + (rand(size(xdotsim))-.5)*noise*max(xdotsim); 
    xdotdotsim = xdotdotsim + (rand(size(xdotdotsim))-.5)*noise*max(xdotdotsim); 
    store_x2 = [store_x2; x(length(x)/2:end)];
    store_xdot2 = [store_xdot2; xdotsim(length(x)/2:end)];
    store_xdotdot2 = [store_xdotdot2; xdotdotsim(length(x)/2:end)];
end

figure(1)

data = [store_x1,store_xdot1;store_x2,store_xdot2];

% call the routine
[bandwidth,density,X,Y]=kde2d(data);
% plot the data and the density estimate
fprintf('bandwidth is %f \n',bandwidth);

figure(1);
%contour3(X,Y,density,50)
hold on
g = scatter(store_x1,store_xdot1,'rx');
hold on

f = plot(store_x2,store_xdot2,'gx');

hold off
title('phase for two tissues')
xlabel('x')
ylabel('xdot')
legend('tissue 1', 'tissue 2')