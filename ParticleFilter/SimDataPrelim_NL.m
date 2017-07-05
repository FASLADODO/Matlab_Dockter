%% Non linear sim

clear all

%settings
fontS = 14;
scale_factor = 0.05;
paramOptions = 2;

%run variables
runs = 10;
paramNoise = 0.1;
dataNoise = 0.01;
trainIDX1 = [1:7,9,10];
trainIDX2 = [1:7,9,10];
runIDX = [8];



%start parameters
phiNL1 = [100;20;1];
phiNL2 = [75;30;1];
phiNL3 = [120;25;1];

type = {'Class1','Class2','Class3'};

%inputs
A = 100;
tend = 10.25;
T = 0.001; % sampling period is fronm 1KHz
t = 0:0.001:tend;

%input is linear increase with time
input.time = t;
input.signals.values = A*t;

input.time = [input.time]';
input.signals.values = [input.signals.values]';
input.signals.dimensions = 1;


fprintf('Running NL simulation 1...')
for kk = 1 :runs

    for jj = 1:length(phiNL1)
        phi(jj) = phiNL1(jj) +(rand(1)-.5)*paramNoise*max(phiNL1(jj));
    end

    % SIMULINK

    sim('SimpleModel_NL.slx');

    % no noise, get output stuff
    u1_tmp = input_out.Data(:,1);
    x1_tmp = state.Data(:,3);
    xdot1_tmp = state.Data(:,2);
    xdotdot1_tmp = state.Data(:,1); 
    state1 = state.Data;

    %Store
    simd_NL{1}.input(:,kk) = u1_tmp;
    simd_NL{1}.state_k{1}(:,kk) = x1_tmp;
    simd_NL{1}.state_k{2}(:,kk) = xdot1_tmp;
    simd_NL{1}.state_k{3}(:,kk) = xdotdot1_tmp;
    simd_NL{1}.params(kk,:) = phi;

end

%Merge all runs
simd_NL{1}.state{1} = reshape(simd_NL{1}.state_k{1}(:,trainIDX1),[],1);
simd_NL{1}.state{2} = reshape(simd_NL{1}.state_k{2}(:,trainIDX1),[],1);
simd_NL{1}.state{3} = reshape(simd_NL{1}.state_k{3}(:,trainIDX1),[],1);

fprintf('Running NL simulation 2...')
for kk = 1 :runs

    for jj = 1:length(phiNL2)
        phi(jj) = phiNL2(jj) +(rand(1)-.5)*paramNoise*max(phiNL2(jj));
    end

    % SIMULINK

    sim('SimpleModel_NL.slx');

    % no noise, get output stuff
    u2_tmp = input_out.Data(:,1);
    x2_tmp = state.Data(:,3);
    xdot2_tmp = state.Data(:,2);
    xdotdot2_tmp = state.Data(:,1); 
    state2 = state.Data;

    %Store
    simd_NL{2}.input(:,kk) = u2_tmp;
    simd_NL{2}.state_k{1}(:,kk) = x2_tmp;
    simd_NL{2}.state_k{2}(:,kk) = xdot2_tmp;
    simd_NL{2}.state_k{3}(:,kk) = xdotdot2_tmp;
    simd_NL{2}.params(kk,:) = phi;

end

%Merge all runs
simd_NL{2}.state{1} = reshape(simd_NL{2}.state_k{1}(:,trainIDX2),[],1);
simd_NL{2}.state{2} = reshape(simd_NL{2}.state_k{2}(:,trainIDX2),[],1);
simd_NL{2}.state{3} = reshape(simd_NL{2}.state_k{3}(:,trainIDX2),[],1);

figure(2)
h1 = quiver(simd_NL{1}.state{1},simd_NL{1}.state{2},simd_NL{1}.state{2}*scale_factor,simd_NL{1}.state{3}*scale_factor,'color',[1 0 0],'AutoScale','off');
hold on
h2 = quiver(simd_NL{2}.state{1},simd_NL{2}.state{2},simd_NL{2}.state{2}*scale_factor,simd_NL{2}.state{3}*scale_factor,'color',[0 1 0],'AutoScale','off');
hold on

title('Non linear model phase portrait*, Two Classes','FontSize',fontS)
xlabel('x','FontSize',fontS)
ylabel('x dot','FontSize',fontS)
legend([h1,h2],type{1},type{2},'FontSize',fontS)

online.state{1} = simd_NL{2}.state_k{1}(:,runIDX);
online.state{2} = simd_NL{2}.state_k{2}(:,runIDX);
online.state{3} = simd_NL{2}.state_k{3}(:,runIDX);

% h3 = quiver(online.state{1},online.state{2},online.state{2}*scale_factor,online.state{3}*scale_factor,'color',[0 0 1],'AutoScale','off');
% hold on

%[ weights,regions,x1grid,x2grid] = DPP_Discriminant_Train_V2( simd_NL, 6, [1,2,3], 1 );


%%
[  testClass, arrClass, listClass  ] = DPP_Discriminant_Online( simd_NL, online, regions, weights, x1grid, x2grid, 5 );




