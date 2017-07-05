
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
paramNoise = 0.01;
dataNoise = 0.02;
trainIDX1 = [1:7,9,10];
trainIDX2 = [1:7,9,10];
runIDX = [8];

%|xdot    | = |0     1|*|x   | + U
%|xdotdot |   |a1 a2| |xdot|

M1 = 2;
D1 = 28;
K1 = 45;
M2 = 2;
D2 = 29;
K2 = 47;


A1 = [0, 1;
     -K1/M1, -D1/M1];
 
A2 = [0, 1;
     -K2/M2, -D2/M2];
     

%Different Parameter vectors
phi1 = [A1(2,1);A1(2,2)];
phi2 = [A2(2,1);A2(2,2)];

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
    
    sim('SimpleModel2.slx');

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
simd{1}.state{1} = reshape(simd{1}.state_k{1}(:,trainIDX1),[],1);
simd{1}.state{2} = reshape(simd{1}.state_k{2}(:,trainIDX1),[],1);
simd{1}.state{3} = reshape(simd{1}.state_k{3}(:,trainIDX1),[],1);
simd{1}.input = reshape(simd{1}.input_k(:,trainIDX1),[],1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%now with second phi%%%%%%%%%%%%%%%%%%%%%%
fprintf('Running simulation 2...')
for kk = 1 :runs
    
    for jj = 1:length(phi2)
        phi(jj) = phi2(jj) +(rand(1)-.5)*paramNoise*max(phi2(jj));
    end

    % SIMULINK
    
    sim('SimpleModel2.slx');
    

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
simd{2}.state{1} = reshape(simd{2}.state_k{1}(:,trainIDX2),[],1);
simd{2}.state{2} = reshape(simd{2}.state_k{2}(:,trainIDX2),[],1);
simd{2}.state{3} = reshape(simd{2}.state_k{3}(:,trainIDX2),[],1);
simd{2}.input = reshape(simd{2}.input_k(:,trainIDX2),[],1);

%%


figure(2)
h1 = quiver(simd{1}.state{2},simd{1}.state{3},simd{1}.state{2}*scale_factor,simd{1}.state{3}*scale_factor,'color',[1 0 0],'AutoScale','off');
hold on
h2 = quiver(simd{2}.state{2},simd{2}.state{3},simd{2}.state{2}*scale_factor,simd{2}.state{3}*scale_factor,'color',[0 1 0],'AutoScale','off');
hold on

online{1}.state{1} = reshape(simd{1}.state_k{1}(:,runIDX),[],1);
online{1}.state{2} = reshape(simd{1}.state_k{2}(:,runIDX),[],1);
online{1}.state{3} = reshape(simd{1}.state_k{3}(:,runIDX),[],1);
online{1}.input = reshape(simd{1}.input_k(:,runIDX),[],1);

online{2}.state{1} = reshape(simd{2}.state_k{1}(:,runIDX),[],1);
online{2}.state{2} = reshape(simd{2}.state_k{2}(:,runIDX),[],1);
online{2}.state{3} = reshape(simd{2}.state_k{3}(:,runIDX),[],1);
online{2}.input = reshape(simd{2}.input_k(:,runIDX),[],1);

h3 = quiver(online{1}.state{2},online{1}.state{3},online{1}.state{2}*scale_factor,online{1}.state{3}*scale_factor,'color',[0 0 1],'AutoScale','off');
hold on

title('linear model phase portrait*, Two Classes','FontSize',fontS)
xlabel('x','FontSize',fontS)
ylabel('x dot','FontSize',fontS)
legend([h1,h2,h3],type{1},type{2},'online','FontSize',fontS)

%get params estimate ----- returns [K;D;M] !!!!!!!!!!
fprintf('TLS estimate 1: \n')
params1_L = TotalLeastSquares( simd,1 )
fprintf('\n Compared to  \n')
mean(simd{1}.params(trainIDX1,:) )

fprintf(' \n TLS estimate 2: \n')
params2_L = TotalLeastSquares( simd,2 )
fprintf('\n Compared to  \n')
mean(simd{2}.params(trainIDX2,:) )


%%

lambdas = [0:0.1:1];

% ----- returns [K;D;M]
cnt = 1;
for L = lambdas
    [ param1(:,cnt), param2(:,cnt) ] = DLS_Discriminant_Train( simd, L );
    lambdars(cnt) = L;
    cnt = cnt + 1;
end

figure(4)
plot(param1(1,:),param1(2,:),'gx','MarkerSize',10)
hold on
plot(param2(1,:),param2(2,:),'kx','MarkerSize',10)
hold on
plot(params1_L(1),params1_L(2),'bo','MarkerSize',10)
hold on
plot(params2_L(1),params2_L(2),'co','MarkerSize',10)
hold off

xlabel('K/M')
ylabel('D/M')
title('Discriminant vs Actual parameters')
legend('param1 d','param2 d','params 1 actual','params 2 actual')

LL = lambdas;

figure(5)
plot(LL,param1(1,:),'gx','MarkerSize',10)
hold on
plot(LL,param1(2,:),'rx','MarkerSize',10)
hold on
plot(LL,param2(1,:),'co','MarkerSize',10)
hold on
plot(LL,param2(2,:),'bo','MarkerSize',10)
hold off

xlabel('index')
ylabel('params')
title('Discriminant vs Actual parameters')
legend('param1 1','param1 2','param2 1','param2 2')

%%
figure(7)
h1 = quiver(simd{1}.state{1},simd{1}.state{2},simd{1}.state{2}*scale_factor,simd{1}.state{3}*scale_factor,'color',[1 0 0],'AutoScale','off');
hold on
h2 = quiver(simd{2}.state{1},simd{2}.state{2},simd{2}.state{2}*scale_factor,simd{2}.state{3}*scale_factor,'color',[0 1 0],'AutoScale','off');
hold on

for cc = 1:length(lambdas)-1
    phi = [-abs(param1(1,cc)), -abs(param1(2,cc)) ];
    sim('SimpleModel2.slx');
        % no noise, get output stuff
    u1_tmp = input_out.Data(:,1);
    x1_tmp = state.Data(:,3);
    xdot1_tmp = state.Data(:,2);
    xdotdot1_tmp = state.Data(:,1); 
    state1 = state.Data;
    stash{1}.state{cc} = [x1_tmp,xdot1_tmp,xdotdot1_tmp];
    
    phi = [-abs(param2(1,cc)), -abs(param2(2,cc)) ];
    sim('SimpleModel2.slx');
    % no noise, get output stuff
    u2_tmp = input_out.Data(:,1); 
    x2_tmp = state.Data(:,3); 
    xdot2_tmp = state.Data(:,2); 
    xdotdot2_tmp = state.Data(:,1); 
    state2 = state.Data;
    
    stash{2}.state{cc} = [x2_tmp,xdot2_tmp,xdotdot2_tmp];
    
    if(max(max(stash{1}.state{cc})) < 50000 && max(max(stash{2}.state{cc})) < 50000 )
        h3 = quiver(x1_tmp,xdot1_tmp,xdot1_tmp*scale_factor,xdotdot1_tmp*scale_factor,'color',[cc/length(param1) 0 0],'AutoScale','off');
        hold on
        h4 = quiver(x2_tmp,xdot2_tmp,xdot2_tmp*scale_factor,xdotdot2_tmp*scale_factor,'color',[0 cc/length(param1) 0],'AutoScale','off');
        hold on
    else
        fprintf('lambda = %f blows up \n',lambdas(cc));
    end
end

hold off

title('linear model phase portrait*, Two Classes','FontSize',fontS)
xlabel('x','FontSize',fontS)
ylabel('x dot','FontSize',fontS)
legend([h1,h2],'class 1','class 2')

%%
%[ weights,regions,x1grid,x2grid] = DPP_Discriminant_Train_V2( simd, 7, [1,2,3], 1 );
% 
% [ testC ] = DPP_Discriminant_Online( simd, online, regions, weights, x1grid, x2grid, 5 );
% 



