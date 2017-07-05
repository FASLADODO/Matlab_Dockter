%Linear model sim

clear all
close all
clc

%settings
fontS = 14;
scale_factor = 0.001;
paramOptions = 2;

%run variables
runs = 10;
paramNoise = 0.05;
dataNoise = 0.1;

trainIDX1 = [1:7,9,10];
trainIDX2 = [1:7,9,10];
trainIDX3 = [1:7,9,10];
runIDX = [8];

%|xdot    | = |0     1|*|x   | + U
%|xdotdot |   |a1 a2| |xdot|

     
%Different Parameter vectors
%start parameters
phi1 = [5;20;1];
phi2 = [7;30;1];
phi3 = [12;25;1];

type = {'Class1','Class2','Class3'};

%inputs
A = 100;
tend = 5;
T = 0.01; % sampling period is fronm 1KHz
t = 0:T:tend;

%input is linear increase with time
input.time = t;
input.signals.values = 10*ones(1,length(t)); %A*t; % rand(1,length(t)) - 0.5;
 
input.time = [input.time]';
input.signals.values = [input.signals.values]';
input.signals.dimensions = 1;

%%%%%%%%%%%%%%%%%%%%%%start with first phi%%%%%%%%%%%%%%%%%%%%%%
fprintf('Running simulation 1...')
for kk = 1 :runs
    
    for jj = 1:length(phi1)
        phi(jj) = phi1(jj) +(rand(1)-.5)*paramNoise*max(phi1(jj));
    end

    % SIMULINK
    
    sim('SimpleModel2_NL.slx');

    % no noise, get output stuff
    u1_tmp = input_out.Data(:,1);
    x1_tmp = state.Data(:,3);
    xdot1_tmp = state.Data(:,2);
    xdotdot1_tmp = state.Data(:,1); 
    state1 = state.Data;

    %add in some noise
    x1_tmp = x1_tmp;% + (rand(size(x1_tmp))-.5)*dataNoise; 
    xdot1_tmp = xdot1_tmp;% + (rand(size(xdot1_tmp))-.5)*dataNoise;
    xdotdot1_tmp = xdotdot1_tmp;% + (rand(size(xdotdot1_tmp))-.5)*dataNoise;
    
    %Store
    simd{1}.input_k(:,kk) = u1_tmp;
    simd{1}.state_k{1}(:,kk) = x1_tmp;
    simd{1}.state_k{2}(:,kk) = xdot1_tmp;
    simd{1}.state_k{3}(:,kk) = xdotdot1_tmp;
    simd{1}.params(kk,:) = phi;

end
fprintf(' DONE\n')

%Merge all runs
simd{1}.state(:,1) = reshape(simd{1}.state_k{1}(:,trainIDX1),[],1);
simd{1}.state(:,2) = reshape(simd{1}.state_k{2}(:,trainIDX1),[],1);
simd{1}.state(:,3) = reshape(simd{1}.state_k{3}(:,trainIDX1),[],1);
simd{1}.state(:,4) = reshape(simd{1}.state_k{3}(:,trainIDX1),[],1) .^2;
simd{1}.state(:,5) = reshape(simd{1}.state_k{3}(:,trainIDX1),[],1) .^3;
simd{1}.input = reshape(simd{1}.input_k(:,trainIDX1),[],1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%now with second phi%%%%%%%%%%%%%%%%%%%%%%
fprintf('Running simulation 2...')
for kk = 1 :runs
    
    for jj = 1:length(phi2)
        phi(jj) = phi2(jj) +(rand(1)-.5)*paramNoise*max(phi2(jj));
    end

    % SIMULINK
    
    sim('SimpleModel2_NL.slx');
    

    % no noise, get output stuff
    u2_tmp = input_out.Data(:,1); 
    x2_tmp = state.Data(:,3); 
    xdot2_tmp = state.Data(:,2); 
    xdotdot2_tmp = state.Data(:,1); 
    state2 = state.Data;

    %add in some noise
    x2_tmp = x2_tmp;% + (rand(size(x2_tmp))-.5)*dataNoise; 
    xdot2_tmp = xdot2_tmp;% + (rand(size(xdot2_tmp))-.5)*dataNoise;
    xdotdot2_tmp = xdotdot2_tmp;% + (rand(size(xdotdot2_tmp))-.5)*dataNoise;

    %Store
    simd{2}.input_k(:,kk) = u2_tmp;
    simd{2}.state_k{1}(:,kk) = x2_tmp;
    simd{2}.state_k{2}(:,kk) = xdot2_tmp;
    simd{2}.state_k{3}(:,kk) = xdotdot2_tmp;
    simd{2}.params(kk,:) = phi;
end
fprintf(' DONE\n')

%Merge all runs
simd{2}.state(:,1) = reshape(simd{2}.state_k{1}(:,trainIDX2),[],1);
simd{2}.state(:,2) = reshape(simd{2}.state_k{2}(:,trainIDX2),[],1);
simd{2}.state(:,3) = reshape(simd{2}.state_k{3}(:,trainIDX2),[],1);
simd{2}.state(:,4) = reshape(simd{2}.state_k{3}(:,trainIDX2),[],1) .^2;
simd{2}.state(:,5) = reshape(simd{2}.state_k{3}(:,trainIDX2),[],1) .^3;
simd{2}.input = reshape(simd{2}.input_k(:,trainIDX2),[],1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%now with second phi%%%%%%%%%%%%%%%%%%%%%%
fprintf('Running simulation 3...')
for kk = 1 :runs
    
    for jj = 1:length(phi3)
        phi(jj) = phi3(jj) +(rand(1)-.5)*paramNoise*max(phi3(jj));
    end

    % SIMULINK
    
    sim('SimpleModel2_NL.slx');
    

    % no noise, get output stuff
    u3_tmp = input_out.Data(:,1); 
    x3_tmp = state.Data(:,3); 
    xdot3_tmp = state.Data(:,2); 
    xdotdot3_tmp = state.Data(:,1); 
    state3 = state.Data;

    %add in some noise
    x3_tmp = x3_tmp;% + (rand(size(x2_tmp))-.5)*dataNoise; 
    xdot3_tmp = xdot3_tmp;% + (rand(size(xdot2_tmp))-.5)*dataNoise;
    xdotdot3_tmp = xdotdot3_tmp;% + (rand(size(xdotdot2_tmp))-.5)*dataNoise;

    %Store
    simd{3}.input_k(:,kk) = u3_tmp;
    simd{3}.state_k{1}(:,kk) = x3_tmp;
    simd{3}.state_k{2}(:,kk) = xdot3_tmp;
    simd{3}.state_k{3}(:,kk) = xdotdot3_tmp;
    simd{3}.params(kk,:) = phi;
end
fprintf(' DONE\n')

%Merge all runs
simd{3}.state(:,1) = reshape(simd{3}.state_k{1}(:,trainIDX3),[],1);
simd{3}.state(:,2) = reshape(simd{3}.state_k{2}(:,trainIDX3),[],1);
simd{3}.state(:,3) = reshape(simd{3}.state_k{3}(:,trainIDX3),[],1);
simd{3}.state(:,4) = reshape(simd{3}.state_k{3}(:,trainIDX3),[],1) .^2;
simd{3}.state(:,5) = reshape(simd{3}.state_k{3}(:,trainIDX3),[],1) .^3;
simd{3}.input = reshape(simd{3}.input_k(:,trainIDX3),[],1);

%%
Parameters = DLS_TrainGeneral(simd, [0,0,0])


%%
figure(2)
h1 = quiver(simd{1}.state(:,1),simd{1}.state(:,2),simd{1}.state(:,2)*scale_factor,simd{1}.state(:,3)*scale_factor,'color',[1 0 0],'AutoScale','off');
hold on
h2 = quiver(simd{2}.state(:,1),simd{2}.state(:,2),simd{2}.state(:,2)*scale_factor,simd{2}.state(:,3)*scale_factor,'color',[0 1 0],'AutoScale','off');
hold on
h3 = quiver(simd{3}.state(:,1),simd{3}.state(:,2),simd{3}.state(:,2)*scale_factor,simd{3}.state(:,3)*scale_factor,'color',[0 0 1],'AutoScale','off');
hold off

% online{1}.state{1} = reshape(simd{1}.state_k{1}(:,runIDX),[],1);
% online{1}.state{2} = reshape(simd{1}.state_k{2}(:,runIDX),[],1);
% online{1}.state{3} = reshape(simd{1}.state_k{3}(:,runIDX),[],1);
% online{1}.input = reshape(simd{1}.input_k(:,runIDX),[],1);
% 
% online{2}.state{1} = reshape(simd{2}.state_k{1}(:,runIDX),[],1);
% online{2}.state{2} = reshape(simd{2}.state_k{2}(:,runIDX),[],1);
% online{2}.state{3} = reshape(simd{2}.state_k{3}(:,runIDX),[],1);
% online{2}.input = reshape(simd{2}.input_k(:,runIDX),[],1);
% 
% h3 = quiver(online{1}.state{2},online{1}.state{3},online{1}.state{2}*scale_factor,online{1}.state{3}*scale_factor,'color',[0 0 1],'AutoScale','off');
% hold on

title('linear model phase portrait*, Two Classes','FontSize',fontS)
xlabel('x','FontSize',fontS)
ylabel('x dot','FontSize',fontS)
%legend([h1,h2,h3],type{1},type{2},type{3},'FontSize',fontS)

%% Test New Params (Parameters[]) DOTHIS

for kk = 1:2
    
    phi(kk) = Parameters(:,kk)';

    % SIMULINK
    
    sim('SimpleModel3_NL.slx');

    % no noise, get output stuff
    u1_tmp = input_out.Data(:,1);
    x1_tmp = state.Data(:,3);
    xdot1_tmp = state.Data(:,2);
    xdotdot1_tmp = state.Data(:,1); 
    state1 = state.Data;

    %add in some noise
    x1_tmp = x1_tmp;% + (rand(size(x1_tmp))-.5)*dataNoise; 
    xdot1_tmp = xdot1_tmp;% + (rand(size(xdot1_tmp))-.5)*dataNoise;
    xdotdot1_tmp = xdotdot1_tmp;% + (rand(size(xdotdot1_tmp))-.5)*dataNoise;
    
    %Store
    simd{1}.input_k(:,kk) = u1_tmp;
    simd{1}.state_k{1}(:,kk) = x1_tmp;
    simd{1}.state_k{2}(:,kk) = xdot1_tmp;
    simd{1}.state_k{3}(:,kk) = xdotdot1_tmp;
    simd{1}.params(kk,:) = phi;