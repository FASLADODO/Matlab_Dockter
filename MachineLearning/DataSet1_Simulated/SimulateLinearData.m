%Create Linear Data Set

fsize = 14;

%data points
nn = 1000;
classes = 2;

% Difference in Parameters (kappa)
AlphaDiff = 1.2; %1.2, 1.1

% Parameters
Alpha1 = [5;2];
AlphaAll = [];
AlphaAll{1} = Alpha1;
rng('default')
for gg = 2:classes
    AlphaAll{gg} = AlphaAll{1}*AlphaDiff;
end

%model noise (kappa)
mdl_noise = 0.001; %noise of parameters

%process noise (lambda)
noise_proc = 0.01; %pos, vel, acc

%Observation noise (epsilon)
obs_noise = 0.01;

%number of runs per class
runs = 100;

%% Actually Run Simulation

%storage
AllData = [];
AllLabels = [];

%Segments
SegmentData = [];
SegmentData.key = {'x','xdot','xdotdot','u'};

for ii = 1:runs%loop through each run
    for gg = 1:classes %loop through each class
        
        %amplitude
        A = 2;
        
        % Input to the system, u = force
        tend = 10.0;
        t = linspace(1,10,tend);
        
        input.time = t;
        input.signals.values = A*t; %ramp input

        input.time = [input.time]';
        input.signals.values = [input.signals.values]';
        input.signals.dimensions = 1;
        
        %current parameters
        alphatrue = AlphaAll{gg};
        alpha = alphatrue + randn(2,1)*mdl_noise;

        % SIMULINK
        fprintf('Simulation, run %d, class %d',ii,gg);
        sim('SimpleLinearModel.slx');
        fprintf(' DONE\n')

        % no noise
        usim = input_out.Data(:,1); %force.Data ;
        xsim = state.Data(:,3); %state.Data(:,3) ; 
        xdotsim = state.Data(:,2); %state.Data(:,3) ; 
        xdotdotsim = state.Data(:,1); %state.Data(:,3) ; 
        statesim = state.Data;
        
        %size of the data
        [ns,~] = size(xsim);

        %Add noise to all
        unoise = usim + (rand(ns,1)-0.5)*range(usim)*obs_noise;
        xnoise  = xsim + (rand(ns,1)-0.5)*range(usim)*noise_proc;
        xdotnoise  = xdotsim + (rand(ns,1)-0.5)*range(usim)*noise_proc;
        xdotdotnoise = xdotdotsim + (rand(ns,1)-0.5)*range(usim)*noise_proc;

        %store
        AllData = [AllData; xnoise, xdotnoise, xdotdotnoise, unoise];
        AllLabels = [AllLabels; ones(ns,1)*gg];
        
        %segment
        SegmentData.Class{gg}.Iteration{ii} = [xnoise, xdotnoise, xdotdotnoise, unoise];
        
    end
end

%% Save it

save('LinearData12.mat','SegmentData');

%% Plot

figure
gscatter(AllData(:,1), AllData(:,4),AllLabels)
xlabel('x','FontSize',fsize)
ylabel('u','FontSize',fsize)
[lh,ic,ip,it]=legend('show');
lh.FontSize = fsize;
lh.Location = 'NorthWest';
