%Linear model sim

clear all
close all
clc

%% sim 2D data

%plot settings
fontS = 14;
scale_factor = 0.001;
paramOptions = 2;

ploton = 0;
    
%run variables
runs = 25;
paramNoise = 0.001;
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


%%%%%%%%%%%%%%%%%%%%%%get random parameters [m, b] (2D data
param_bar = [28,45];

%gamma = 0.05;
%param_arr(1,:) = Create_Gamma_Params(param_bar, gamma); 

gamma_arr = [0, 0.01, 0.02, 0.05, 0.1, 0.2];        

for kk = 1:length(gamma_arr)
   param_arr(kk,:) = Create_Gamma_Params1D(param_bar, gamma_arr(kk)); 
   gamma(kk) =  Compute_Gamma(  param_arr(kk,2), param_bar(2) );
end



%run through all runs with comparison parameters (phi1)
A1 = [1/param_bar(1), -param_bar(2)/param_bar(1)];
phi1 = [A1(1),A1(2)];

pp = 1; %just for class 1
fprintf('Running simulation 1...')
for kk = 1 :runs

    for jj = 1:length(phi1)
        phi(jj) = phi1(jj) +(randn(1)-.5)*paramNoise*phi1(jj);
    end

    % SIMULINK

    sim('SimpleModel1D.slx');

    % no noise, get output stuff
    u1_tmp = input_out.Data(:,1);
    x1_tmp = state.Data(:,2);
    xdot1_tmp = state.Data(:,1);
    state1 = state.Data;

    %add in some noise
    x1_tmp = x1_tmp;% + (randn(size(x1_tmp))-.5)*dataNoise*mean(x1_tmp); 
    xdot1_tmp = xdot1_tmp;% + (randn(size(xdot1_tmp))-.5)*dataNoise*mean(xdot1_tmp);
    
    %get errors
    D1 = [xdot1_tmp,x1_tmp];
    param_test = abs(phi)'; %[param_bar(3)/param_bar(1);param_bar(2)/param_bar(1);1];
    e_1_tmp(:,kk) = u1_tmp - (D1*param_test); 
    
%     figure(400+kk)
%     plot(1:length(u1_tmp),u1_tmp,'rx');
%     hold on
%     plot(1:length(D1*param_test),D1*param_test,'bx');
%     hold off

    %Store [xdot,x]
    simd{1}.input_k(:,kk) = u1_tmp;
    simd{1}.state_k{1}(:,kk) = xdot1_tmp;
    simd{1}.state_k{2}(:,kk) = x1_tmp;
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
    avg_x = mean(simd{1}.state_k{2},2);
    avg_xdot = mean(simd{1}.state_k{1},2);
    
    phi = -param_bar(2)/param_bar(1);
    sim('SimpleModel1D.slx');
    
        % no noise, get output stuff
    u1_tmp = input_out.Data(:,1);
    x1_tmp = state.Data(:,2);
    xdot1_tmp = state.Data(:,1);
    state1 = state.Data;
    
    figure(101+pp)
    for kk = 1 :runs
        hh1 = plot(simd{1}.state_k{2}(:,kk),simd{1}.state_k{1}(:,kk),'bx');
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
simd{1}.input = reshape(simd{1}.input_k,[],1);

%Loop through all possible param variations
for pp = 1:length(param_arr)
    
    A2 = [1/param_arr(pp,1),-param_arr(pp,2)/param_arr(pp,1)];

    %Variation Parameter vector
    phi2 = [A2(1),A2(2)];
    
    fprintf('Running simulation 2...')
    for kk = 1 :runs

        for jj = 1:length(phi2)
            phi(jj) = phi2(jj) +(randn(1)-.5)*paramNoise*phi2(jj);
        end

        % SIMULINK

        sim('SimpleModel1D.slx');


        % no noise, get output stuff
        u2_tmp = input_out.Data(:,1); 
        x2_tmp = state.Data(:,2); 
        xdot2_tmp = state.Data(:,1);  
        state2 = state.Data;

        %add in some noise
        x2_tmp = x2_tmp;% + (randn(size(x2_tmp))-.5)*dataNoise*mean(x2_tmp); 
        xdot2_tmp = xdot2_tmp;% + (randn(size(xdot2_tmp))-.5)*dataNoise*mean(xdot2_tmp);

        %get errors
        D2 = [xdot2_tmp,x2_tmp];
        param_test = abs(phi)'; %[param_arr(pp,3)/param_arr(pp,1);param_arr(pp,2)/param_arr(pp,1);1];
        e_2_tmp(:,kk) = u2_tmp - (D2*param_test); 

        
        %Store [xdot,x]
        simd{2}.gamma{pp}.input_k(:,kk) = u2_tmp;
        simd{2}.gamma{pp}.state_k{1}(:,kk) = xdot2_tmp;
        simd{2}.gamma{pp}.state_k{2}(:,kk) = x2_tmp;
        simd{2}.gamma{pp}.params(kk,:) = phi;
    end
    fprintf(' DONE\n')
   
    %compute class inherint noise
    ec_2(:,pp) = mean(abs(e_2_tmp),2) ./2;
    
    %Merge all runs
    simd{2}.gamma{pp}.state{1} = reshape(simd{2}.gamma{pp}.state_k{1},[],1);
    simd{2}.gamma{pp}.state{2} = reshape(simd{2}.gamma{pp}.state_k{2},[],1);
    simd{2}.gamma{pp}.input = reshape(simd{2}.gamma{pp}.input_k,[],1);
    
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


h1 = scatter(simd{1}.state{2},simd{1}.state{1},[],[0 0 0]);
hold on

pp = 1;
shape1 = reshape(simd{2}.gamma{pp}.state_k{1},[],1);
datar1 = shape1(1:stepps:end,:);
shape2 = reshape(simd{2}.gamma{pp}.state_k{2},[],1);
datar2 = shape2(1:stepps:end,:);
h(pp) = scatter(datar2,datar1,'cx','LineWidth',0.2);
hold on 
pp = 2;
shape1 = reshape(simd{2}.gamma{pp}.state_k{1},[],1);
datar1 = shape1(1:stepps:end,:);
shape2 = reshape(simd{2}.gamma{pp}.state_k{2},[],1);
datar2 = shape2(1:stepps:end,:);
h(pp) = scatter(datar2,datar1,'rx','LineWidth',0.5);
hold on 
pp = 3;
shape1 = reshape(simd{2}.gamma{pp}.state_k{1},[],1);
datar1 = shape1(1:stepps:end,:);
shape2 = reshape(simd{2}.gamma{pp}.state_k{2},[],1);
datar2 = shape2(1:stepps:end,:);
h(pp) = scatter(datar2,datar1,'gx','LineWidth',0.5);
hold off



str_p=sprintf('Linear model phase portraits');

title(str_p,'FontSize',fontS)
xlabel('x','FontSize',fontS)
ylabel('x dot','FontSize',fontS)
h_legend1=legend([h1,h(1),h(2),h(3)],'class 1','class 2 \gamma=0','class 2 \gamma=0.01','class 2 \gamma=0.02');
set(h_legend1,'FontSize',12);


%[0, 0.01, 0.02, 0.05, 0.1, 0.2];

%% Get all relevant states and lambdas

%Compute numeric version of xbar and ubar and lambda critical
for ii = 1:length(param_arr)
    datatemp{1} = simd{1};
    datatemp{2} = simd{2}.gamma{ii};
    
    %XBar and UBar
    States_All{ii}.Mod_States = ModifyStates(datatemp);
    
    %general params
    States_All{ii}.LambdaCrit = Lambda_Critical(datatemp);
    
    %general params
    States_All{ii}.ParamG = GenericParameters(datatemp);

end

%% Get all the discriminant parameters

percentz = 0:0.01:0.99;


%Compute discriminant parameters
for ii = 1:length(param_arr)
    datatemp{1} = simd{1};
    datatemp{2} = simd{2}.gamma{ii};
    
    for pp = 1:length(percentz)
        DLS_Params{ii,pp} = DLS_TrainGeneral_LambdaCritical(datatemp, States_All{ii}.LambdaCrit, percentz(pp));
        
        tempr1{ii}(:,pp) = DLS_Params{ii,pp}{1,2};
        tempr2{ii}(:,pp) = DLS_Params{ii,pp}{2,1};
    end
    
    
    
    %plot all dem params 
    hFig = figure(ii);
    wid = 640;
    hei = 1000;
    set(hFig, 'Position', [(ii-1)*wid 0 wid hei])

    %class 1
    subplot(4,1,1)
    scatter( percentz, tempr1{ii}(1,:), 'r.');
    hold on
    plot( [percentz(1),percentz(end)], [States_All{ii}.ParamG(1,1),States_All{ii}.ParamG(1,1)], 'b' );
    hold off
    xlabel('\zeta','FontSize',fontS)
    ylabel('Phi 1','FontSize',fontS)
    
    titString = sprintf('Change in Parameters Class: %d, gamma: %f', 1, gamma_arr(ii));
    title(titString)
    
    subplot(4,1,2)
    scatter( percentz, tempr1{ii}(2,:), 'r.');
    hold on
    plot( [percentz(1),percentz(end)], [States_All{ii}.ParamG(2,1),States_All{ii}.ParamG(2,1)], 'b' );
    hold off
    xlabel('\zeta','FontSize',fontS)
    ylabel('Phi 2','FontSize',fontS)

    %class 2
    subplot(4,1,3)
    scatter( percentz, tempr2{ii}(1,:), 'r.');
    hold on
    plot( [percentz(1),percentz(end)], [States_All{ii}.ParamG(1,2),States_All{ii}.ParamG(1,2)], 'b' );
    hold off
    xlabel('\zeta','FontSize',fontS)
    ylabel('Phi 1','FontSize',fontS)
    
    titString = sprintf('Change in Parameters Class: %d, gamma: %f', 2, gamma_arr(ii));
    title(titString)

    subplot(4,1,4)
    scatter( percentz, tempr2{ii}(2,:), 'r.');
    hold on
    plot( [percentz(1),percentz(end)], [States_All{ii}.ParamG(2,2),States_All{ii}.ParamG(2,2)], 'b' );
    hold off
    
    xlabel('\zeta','FontSize',fontS)
    ylabel('Phi 2','FontSize',fontS)
    legend('DLS params','true parameter')
end

%% Leave one out classifcation

classes = 2;

%can vary this
lamdaper = 0.95;

totalacc = zeros(1,length(param_arr));
plotacc1 = [];
plotacc2 = [];

for pp = 1:length(param_arr)
    for kk = 1 :runs
        %Merge all but one runs for training
        temp1 = simd{1}.state_k{1};
        temp2 = simd{1}.state_k{2};
        tempin = simd{1}.input_k;
        temp1(:,kk) = [];
        temp2(:,kk) = [];
        tempin(:,kk) = [];
        dataLOO{1}.state{1} = reshape(temp1,[],1);
        dataLOO{1}.state{2} = reshape(temp2,[],1);
        dataLOO{1}.input = reshape(tempin,[],1);

        temp1 = simd{2}.gamma{pp}.state_k{1};
        temp2 = simd{2}.gamma{pp}.state_k{2};
        tempin = simd{2}.gamma{pp}.input_k;
        temp1(:,kk) = [];
        temp2(:,kk) = [];
        tempin(:,kk) = [];
        dataLOO{2}.state{1} = reshape(temp1,[],1);
        dataLOO{2}.state{2} = reshape(temp2,[],1);
        dataLOO{2}.input = reshape(tempin,[],1);

        %compute LOO params
        for ll = 1:length(percentz)
            DLS_Param_LOO{pp,kk}.percent{ll} = DLS_TrainGeneral_LambdaCritical(dataLOO, States_All{pp}.LambdaCrit, percentz(ll));
        end
        
        %now grab useen data
        dataOnline{1}.state{1} = simd{1}.state_k{1}(:,kk);
        dataOnline{1}.state{2} = simd{1}.state_k{2}(:,kk);
        dataOnline{1}.input = simd{1}.input_k(:,kk);
        
        dataOnline{2}.state{1} = simd{2}.gamma{pp}.state_k{1}(:,kk);
        dataOnline{2}.state{2} = simd{2}.gamma{pp}.state_k{2}(:,kk);
        dataOnline{2}.input = simd{2}.gamma{pp}.input_k(:,kk);
        
        %try classification with unseen data
        for jj=1:classes
            %try all percents for now
            for ll = 1:length(percentz)
                Classify{jj}.gamma{pp}.percent{ll}.run{kk} = DLS_Online_LOO(dataOnline, jj, DLS_Param_LOO{pp,kk}.percent{ll} );

                %get accuracy
                acctemp = Classify{jj}.gamma{pp}.percent{ll}.run{kk}.classes(:,1) == Classify{jj}.gamma{pp}.percent{ll}.run{kk}.classes(:,2);
                Accuracy{jj}.gamma{pp}.percent{ll}.run(kk) = sum(acctemp)/length(acctemp);
            end
        end
    end
    
    %save accuracy for all percents for later
    for ll = 1:length(percentz)
        totalacc(pp,ll) = mean([Accuracy{1}.gamma{pp}.percent{ll}.run,Accuracy{2}.gamma{pp}.percent{ll}.run]);
        %for plots
        plotacc1 = [ plotacc1; gamma_arr(pp),percentz(ll), mean(Accuracy{1}.gamma{pp}.percent{ll}.run)];
        plotacc2 = [ plotacc2; gamma_arr(pp),percentz(ll), mean(Accuracy{2}.gamma{pp}.percent{ll}.run)];
    end
end


%% Plot

fontS = 14;

figure
c1 = Surface3D(plotacc1(:,1),plotacc1(:,2),plotacc1(:,3),'mesh');
xlabel('\gamma','FontSize',fontS)
ylabel('\zeta','FontSize',fontS)
zlabel('Accuracy','FontSize',fontS)
title('Accuracy vs \gamma and \zeta, Class 1')

figure
c2 = Surface3D(plotacc2(:,1),plotacc2(:,2),plotacc2(:,3),'mesh');
xlabel('\gamma','FontSize',fontS)
ylabel('\zeta','FontSize',fontS)
zlabel('Accuracy','FontSize',fontS)
title('Accuracy vs \gamma and \zeta, Class 2')

