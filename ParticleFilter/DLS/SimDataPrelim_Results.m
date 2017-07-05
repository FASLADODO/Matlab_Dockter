
%Linear model sim

clear all
close all
clc

%%

%plot settings
fontS = 14;
scale_factor = 0.001;
paramOptions = 2;

ploton = 0;
    
%run variables
runs = 25;
paramNoise = 0.01;
dataNoise = 0.00001;
trainIDX1 = [1:7,9,10];
trainIDX2 = [1:7,9,10];
runIDX = [8];

%%%%%%%%%%%% inputs
A = 100;
tend = 8;
T = 0.01; % sampling period is fronm 1KHz
t = 0:T:tend;

%input is linear increase with time
input.time = t;
input.signals.values = 10*ones(1,length(t));%A*t; % rand(1,length(t)) - 0.5;

input.time = [input.time]';
input.signals.values = [input.signals.values]';
input.signals.dimensions = 1;


%|xdot    | = |0     1|*|x   | + U
%|xdotdot |   |a1 a2| |xdot|

lambdas0_1 = [0:0.05:1]; %%%%Only in this range
% lambdas1_10 = [1:0.5:10];
% lambdas0_10 = [0:0.1:1,1.5:0.5:10];


%%%%%%%%%%%%%%%%%%%%%%get random parameters [M,D,K]
param_bar = [2,28,45];

% param_arr = [5,12,23;
%             5,14,25;
%             5,20,30;
%             5,25,45;
%             5,50,60;
%             5,100,250];
gamma_arr = [0, 0.01, 0.02, 0.05, 0.1, 0.2];        

for kk = 1:length(gamma_arr)
   param_arr(kk,:) = Create_Gamma_Params(param_bar, gamma_arr(kk)); 
   gamma(kk) =  Compute_Gamma(  param_arr(kk,2:3), param_bar(2:3) );
end


%run through all runs with comparison parameters (phi1)
A1 = [0, 1;
     -param_bar(3)/param_bar(1), -param_bar(2)/param_bar(1)];
phi1 = [A1(2,1);A1(2,2)];

pp = 1; %just for class 1
fprintf('Running simulation 1...')
for kk = 1 :runs

    for jj = 1:length(phi1)
        phi(jj) = phi1(jj) +(randn(1)-.5)*paramNoise*phi1(jj);
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
    x1_tmp = x1_tmp;% + (randn(size(x1_tmp))-.5)*dataNoise*mean(x1_tmp); 
    xdot1_tmp = xdot1_tmp;% + (randn(size(xdot1_tmp))-.5)*dataNoise*mean(xdot1_tmp);
    xdotdot1_tmp = xdotdot1_tmp ;%+ (randn(size(xdotdot1_tmp))-.5)*dataNoise*mean(xdotdot1_tmp);
    
    %get errors
    D1 = [x1_tmp,xdot1_tmp,xdotdot1_tmp];
    param_test = abs([phi(1);phi(2);1]); %[param_bar(3)/param_bar(1);param_bar(2)/param_bar(1);1];
    e_1_tmp(:,kk) = u1_tmp - (D1*param_test); 
    
%     figure(400+kk)
%     plot(1:length(u1_tmp),u1_tmp,'rx');
%     hold on
%     plot(1:length(D1*param_test),D1*param_test,'bx');
%     hold off

    %Store
    simd{1}.input_k(:,kk) = u1_tmp;
    simd{1}.state_k{1}(:,kk) = x1_tmp;
    simd{1}.state_k{2}(:,kk) = xdot1_tmp;
    simd{1}.state_k{3}(:,kk) = xdotdot1_tmp;
    simd{1}.params(kk,:) = phi;

end
fprintf(' DONE\n')

%compute class inherint noise
ec_1(:,pp) = mean(abs(e_1_tmp),2);

if(ploton)
    figure(100+pp)
    plot(1:length(ec_1(:,pp)),ec_1(:,pp),'rx');
    
    title('inherent system noise (class 1)')
    xlabel('time')
    ylabel('error')
    
    %See noises
    avg_x = mean(simd{1}.state_k{1},2);
    avg_xdot = mean(simd{1}.state_k{2},2);
    avg_xdotdot = mean(simd{1}.state_k{3},2);
    
    phi = [-param_bar(3)/param_bar(1), -param_bar(2)/param_bar(1) ];
    sim('SimpleModel2.slx');
        % no noise, get output stuff
    u1_tmp = input_out.Data(:,1);
    x1_tmp = state.Data(:,3);
    xdot1_tmp = state.Data(:,2);
    xdotdot1_tmp = state.Data(:,1); 
    state1 = state.Data;
    
    figure(101+pp)
    for kk = 1 :runs
        hh1 = plot(simd{1}.state_k{1}(:,kk),simd{1}.state_k{2}(:,kk),'bx');
        hold on
    end
    hh2 = plot(x1_tmp,xdot1_tmp,'gx');
    hold on
    hh3 = plot(avg_x+abs(ec_1),avg_xdot+abs(ec_1),'rx');
    hold on
    hh4 = plot(avg_x-abs(ec_1),avg_xdot-abs(ec_1),'rx');
    hold off
    
    title('Simulated data with error bars','FontSize',fontS)
    xlabel('x','FontSize',fontS)
    ylabel('xdot','FontSize',fontS)
    legend([hh1,hh2,hh3],'simdata','Actual','error bars','FontSize',12)
end

%Merge all class 1 runs for comparison
simd{1}.state{1} = reshape(simd{1}.state_k{1},[],1);
simd{1}.state{2} = reshape(simd{1}.state_k{2},[],1);
simd{1}.state{3} = reshape(simd{1}.state_k{3},[],1);
simd{1}.input = reshape(simd{1}.input_k,[],1);

%Loop through all possible param variations
for pp = 1:length(param_arr)

    A2 = [0, 1;
         -param_arr(pp,3)/param_arr(pp,1), -param_arr(pp,2)/param_arr(pp,1)];


    %Variation Parameter vector
    phi2 = [A2(2,1);A2(2,2)];
    
    fprintf('Running simulation 2...')
    for kk = 1 :runs

        for jj = 1:length(phi2)
            phi(jj) = phi2(jj) +(randn(1)-.5)*paramNoise*phi2(jj);
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
        x2_tmp = x2_tmp;% + (randn(size(x2_tmp))-.5)*dataNoise*mean(x2_tmp); 
        xdot2_tmp = xdot2_tmp;% + (randn(size(xdot2_tmp))-.5)*dataNoise*mean(xdot2_tmp);
        xdotdot2_tmp = xdotdot2_tmp;% + (randn(size(xdotdot2_tmp))-.5)*dataNoise;%*max(xdotdot2_tmp);

        %get errors
        D2 = [x2_tmp,xdot2_tmp,xdotdot2_tmp];
        param_test = abs([phi(1);phi(2);1]); %[param_arr(pp,3)/param_arr(pp,1);param_arr(pp,2)/param_arr(pp,1);1];
        e_2_tmp(:,kk) = u2_tmp - (D2*param_test); 

        
        %Store
        simd{2}.gamma{pp}.input_k(:,kk) = u2_tmp;
        simd{2}.gamma{pp}.state_k{1}(:,kk) = x2_tmp;
        simd{2}.gamma{pp}.state_k{2}(:,kk) = xdot2_tmp;
        simd{2}.gamma{pp}.state_k{3}(:,kk) = xdotdot2_tmp;
        simd{2}.gamma{pp}.params(kk,:) = phi;
    end
    fprintf(' DONE\n')
   
    %compute class inherint noise
    ec_2(:,pp) = mean(abs(e_2_tmp),2) ./2;
    
    %Merge all runs
    simd{2}.state{1} = reshape(simd{2}.gamma{pp}.state_k{1},[],1);
    simd{2}.state{2} = reshape(simd{2}.gamma{pp}.state_k{2},[],1);
    simd{2}.state{3} = reshape(simd{2}.gamma{pp}.state_k{3},[],1);
    simd{2}.input = reshape(simd{2}.gamma{pp}.input_k,[],1);

    if(ploton)
        figure(pp)
        h1 = quiver(simd{1}.state{1},simd{1}.state{2},simd{1}.state{2}*scale_factor,simd{1}.state{3}*scale_factor,'color',[1 0 0],'AutoScale','off');
        hold on
        h2 = quiver(simd{2}.state{1},simd{2}.state{2},simd{2}.state{2}*scale_factor,simd{2}.state{3}*scale_factor,'color',[0 1 0],'AutoScale','off');
        hold on

        str_p=sprintf('linear model phase portrait, gamma = %f ', gamma(pp));

        title(str_p,'FontSize',fontS)
        xlabel('x','FontSize',fontS)
        ylabel('x dot','FontSize',fontS)
        legend([h1,h2],'class 1','class 2')
        
        figure(30+pp)
        h1 = quiver(simd{1}.state{2},simd{1}.state{3},simd{1}.state{2}*scale_factor,simd{1}.state{3}*scale_factor,'color',[1 0 0],'AutoScale','off');
        hold on
        h2 = quiver(simd{2}.state{2},simd{2}.state{3},simd{2}.state{2}*scale_factor,simd{2}.state{3}*scale_factor,'color',[0 1 0],'AutoScale','off');
        hold on

        str_p=sprintf('linear model phase portrait, gamma = %f ', gamma(pp));

        title(str_p,'FontSize',fontS)
        xlabel('x dot','FontSize',fontS)
        ylabel('x dotdot','FontSize',fontS)
        legend([h1,h2],'class 1','class 2')
        
        figure(300 + pp)
        plot(1:length(ec_2(:,pp)),ec_2(:,pp),'rx');
        title('inherent system noise  (class 2)')
        xlabel('time')
        ylabel('error')
        
    end
    
end

%% Plot phase portraits for each class

colormap = {[1 0 0],[0 1 0], [0 0 1], [0 1 1], [1 0 1], [0 1/2 1/2] , [1/2 0 1]};

% '--gs',...
%     'LineWidth',2,...
%     'MarkerSize',10,...
%     'MarkerEdgeColor','b',...
%     'MarkerFaceColor',[0.5,0.5,0.5])

figure(111)

stepps = 25;


h1 = quiver(simd{1}.state{1},simd{1}.state{2},simd{1}.state{2}*0.01,simd{1}.state{3}*0.01,'color',[0 0 0],'AutoScale','off');
hold on

pp = 3;
shape1 = reshape(simd{2}.gamma{pp}.state_k{1},[],1);
datar1 = shape1(1:stepps:end,:);
shape2 = reshape(simd{2}.gamma{pp}.state_k{2},[],1);
datar2 = shape2(1:stepps:end,:);
h(pp) = scatter(datar1,datar2,'cx','LineWidth',0.2);
hold on 
pp = 5;
shape1 = reshape(simd{2}.gamma{pp}.state_k{1},[],1);
datar1 = shape1(1:stepps:end,:);
shape2 = reshape(simd{2}.gamma{pp}.state_k{2},[],1);
datar2 = shape2(1:stepps:end,:);
h(pp) = scatter(datar1,datar2,'rx','LineWidth',0.5);
hold on 
pp = 6;
shape1 = reshape(simd{2}.gamma{pp}.state_k{1},[],1);
datar1 = shape1(1:stepps:end,:);
shape2 = reshape(simd{2}.gamma{pp}.state_k{2},[],1);
datar2 = shape2(1:stepps:end,:);
h(pp) = scatter(datar1,datar2,'gx','LineWidth',0.5);
hold off



str_p=sprintf('Linear model phase portraits');

title(str_p,'FontSize',fontS)
xlabel('x','FontSize',fontS)
ylabel('x dot','FontSize',fontS)
h_legend1=legend([h1,h(3),h(5),h(6)],'class 1','class 2 \gamma=0.02','class 2 \gamma=0.1','class 2 \gamma=0.2');
set(h_legend1,'FontSize',12);


%[0, 0.01, 0.02, 0.05, 0.1, 0.2];

%% Test online classification using leave one out

class_actual = 1;

lambdas = lambdas0_1;
%pool=parpool(8);
tic;
parfor pp = 1:length(param_arr)

    classifications{pp} = LeaveOneOut_Par( pp, simd, lambdas, runs, class_actual, ec_1, ec_2 );

end
toc;
%% Plot new trajectories for each parameters as lambda increases

pp = 6;

simd{2}.state{1} = reshape(simd{2}.gamma{pp}.state_k{1},[],1);
simd{2}.state{2} = reshape(simd{2}.gamma{pp}.state_k{2},[],1);
simd{2}.state{3} = reshape(simd{2}.gamma{pp}.state_k{3},[],1);
simd{2}.input = reshape(simd{2}.gamma{pp}.input_k,[],1);

figure(pp)
h1 = quiver(simd{1}.state{1},simd{1}.state{2},simd{1}.state{2}*scale_factor,simd{1}.state{3}*scale_factor,'color',[1 0 1],'AutoScale','off');
hold on
h2 = quiver(simd{2}.state{1},simd{2}.state{2},simd{2}.state{2}*scale_factor,simd{2}.state{3}*scale_factor,'color',[0 1 1],'AutoScale','off');
hold on

rr = 1;
for cc = 1:length(lambdas)-1
    %phi = [-abs(classifications{pp}.lambda{cc}.params(1,rr)), -abs(classifications{pp}.lambda{cc}.params(2,rr)) ];
    phi = [-classifications{pp}.lambda{cc}.params(1,rr), -classifications{pp}.lambda{cc}.params(2,rr) ];
    
    sim('SimpleModel2.slx');
        % no noise, get output stuff
    u1_tmp = input_out.Data(:,1);
    x1_tmp = state.Data(:,3);
    xdot1_tmp = state.Data(:,2);
    xdotdot1_tmp = state.Data(:,1); 
    state1 = state.Data;
    stash{1}.state{cc} = [x1_tmp,xdot1_tmp,xdotdot1_tmp];
    
    %phi = [-abs(classifications{pp}.lambda{cc}.params(4,rr)), -abs(classifications{pp}.lambda{cc}.params(5,rr)) ];
    phi = [-classifications{pp}.lambda{cc}.params(4,rr), -classifications{pp}.lambda{cc}.params(5,rr) ];
    
    sim('SimpleModel2.slx');
    % no noise, get output stuff
    u2_tmp = input_out.Data(:,1); 
    x2_tmp = state.Data(:,3); 
    xdot2_tmp = state.Data(:,2); 
    xdotdot2_tmp = state.Data(:,1); 
    state2 = state.Data;
    
    stash{2}.state{cc} = [x2_tmp,xdot2_tmp,xdotdot2_tmp];
    
    if(max(max(stash{1}.state{cc})) < 50000 && max(max(stash{2}.state{cc})) < 50000 )
        h3 = quiver(x1_tmp,xdot1_tmp,xdot1_tmp*scale_factor,xdotdot1_tmp*scale_factor,'color',[cc/length(lambdas) 0 0],'AutoScale','off');
        hold on
        h4 = quiver(x2_tmp,xdot2_tmp,xdot2_tmp*scale_factor,xdotdot2_tmp*scale_factor,'color',[0 cc/length(lambdas) 0],'AutoScale','off');
        hold on
    else
        fprintf('lambda = %f blows up \n',lambdas(cc));
    end
end

hold off

title('linear model phase portrait w/ discriminant parameters','FontSize',fontS)
xlabel('x','FontSize',fontS)
ylabel('x dot','FontSize',fontS)
legend([h1,h2,h3,h4],'class 1','class 2','DLS 1','DLS 2')

%% Plot accuracy versus lambda

%load classifcations.mat

colormap = {[1 0 0],[0 1 0], [0 0 1], [0 1 1], [1 0 1], [0 1/2 1] , [1/2 0 1]};

figure(24)
shift = 0.01;
for pp = 1:length(gamma)
    classifications{pp}.percent = [];
    for L = 1:length(lambdas)
        idx_corr = find(classifications{pp}.lambda{L}.class(1,:) == classifications{pp}.lambda{L}.class(2,:));
        percent_correct = length(idx_corr)/length(classifications{pp}.lambda{L}.class(1,:))*100;
        classifications{pp}.percent(L) = percent_correct;
    end
    %figure(pp)
    h(pp) = plot(lambdas+shift*pp,classifications{pp}.percent+shift*pp,'-x','Color', colormap{pp},'LineWidth',3);
    hold on
    str{pp} = sprintf('gamma = %f',gamma(pp));
    %title(str{pp})
    %axis([-5 25 -5 110])
end
hold off

axis([min(lambdas) max(lambdas) -1 101])
title('Discriminant Dynamics Improved Classification','FontSize',fontS)
xlabel('% \lambda_{C}','FontSize',fontS)
ylabel('% Correct Classification','FontSize',fontS)
h_legend=legend(h,str);
set(h_legend,'FontSize',12);







