%make up data using Donor.Tissue.Location.Grasp

%settings
fontS = 14;

%run variables
patientNoise = [0.05,0.05,0.05];
paramNoise = [0.01,0.01,0.01];
dataNoise = [0.002,0.002,0.002];
inputNoise = 0.03;

grasps = 20;
patients = 5;
locations = 3;

%ACTUAL PARAMETERS
%F = D*xdot + a*(exp(-B*x) - 1)
d1 = 0.02;
a1 = 0.1;
b1 = 0.01;

d2 = 0.9;
a2 = 0.1;
b2 = 0.4;


F = 4;
     
%Parameter vectors (d,alpha,beta)
paramsTrue{1} = [d1,a1,b1];
paramsTrue{2} = [d2,a2,b2];

%for classes
key.cslist = {'Class1','Class2'};
key.paramsize = length(paramsTrue{1});
key.theta = 3;
key.thetadot = 2;
key.thetadotdot = 1;
key.dataColumns = [1,2,3];
key.inputColumns = 4;

%storage
GraspData = [];
PatientData = [];
AllData = [];
AllLabels = [];
ClassData{1} = [];
ClassData{2} = [];
ClassInput{1} = [];
ClassInput{2} = [];


%simulink things
tend = 1;
T = 0.001; % sampling period is fronm 1KHz
t = 0:T:tend;

%input is linear increase with time
input.time = t;
input.signals.values = F*t; %10*ones(1,length(t));
input.time = [input.time]';
input.signals.values = [input.signals.values]';
input.signals.dimensions = 1;


%% Run simulink model and store


fprintf('Running simulations  \n')
for pp = 1:patients
    pp
    for gg = 1:length(key.cslist)
        paramsPatient{gg} = paramsTrue{gg} + randn(1,key.paramsize).*patientNoise;
        for ll = 1:locations
            %slight variation for each location
            paramBase = paramsPatient{gg} + randn(1,key.paramsize).*0.0001;
            for ss = 1:grasps
                %these are what actually get used in simulink (true params with a
                %bit of noise)
                phi = paramBase + randn(1,key.paramsize).*paramNoise;

                % SIMULINK
                sim('TissueModel1.slx');

                % no noise, get output stuff
                utrue = input_out.Data(:,1);
                statetrue = state.Data;

                %add in some noise
                unoise = utrue + inputNoise;
                statenoise = statetrue + bsxfun(@times,randn(size(statetrue)),dataNoise);

                %Store individual grasps for all patients
                GraspData{gg}.grasp{ss}.InputTrue = utrue;
                GraspData{gg}.grasp{ss}.StateTrue = statetrue;
                GraspData{gg}.grasp{ss}.Input = unoise;
                GraspData{gg}.grasp{ss}.State = statenoise;
                GraspData{gg}.grasp{ss}.Time = t;
                GraspData{gg}.grasp{ss}.Phi = phi;

                %Store individual grasps for each patient
                SegData.Donor{pp}.Tissue{gg}.Location{ll}.Grasp{ss}.InputTrue = utrue;
                SegData.Donor{pp}.Tissue{gg}.Location{ll}.Grasp{ss}.StateTrue = statetrue;
                SegData.Donor{pp}.Tissue{gg}.Location{ll}.Grasp{ss}.Input = unoise;
                SegData.Donor{pp}.Tissue{gg}.Location{ll}.Grasp{ss}.State = statenoise;
                SegData.Donor{pp}.Tissue{gg}.Location{ll}.Grasp{ss}.Time = t;
                SegData.Donor{pp}.Tissue{gg}.Location{ll}.Grasp{ss}.Phi = phi;

                %class data for stuff
                ClassData{gg} = [ClassData{gg}; statenoise];
                ClassInput{gg} = [ClassInput{gg}; utrue];

                %store all for test plots
                AllData = [AllData; statenoise, unoise];
                AllLabels = [AllLabels; ones(length(unoise),1)*gg];
            end
        end
    end
end
fprintf(' DONE\n')

figure
gscatter(AllData(:,key.theta),AllData(:,key.thetadot),AllLabels)
xlabel('theta')
ylabel('thetadot')

figure
gscatter(AllData(:,key.theta),AllData(:,end),AllLabels)
xlabel('theta')
ylabel('force')

figure
gscatter3(AllData(:,key.theta),AllData(:,key.thetadot),AllData(:,end),AllLabels)
xlabel('theta')
ylabel('thetadot')
zlabel('force')



%% Save this all to a mat

save segData.mat AllData AllLabels ClassData ClassInput GraspData SegData key input

