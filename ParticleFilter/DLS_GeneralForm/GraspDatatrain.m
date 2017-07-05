clear all

load('HardTissue.mat');% Data_Hard 
load('SoftTissue.mat');% Data_Soft

%Area of the load cell we think maybe
AREAZ = 1.86325e-04;

hard_id = 1;
soft_id = 2;


%count number of grasps
num_Hard_Grasps = length(Data_Hard);
num_Soft_Grasps = length(Data_Soft);

ploton = false;

if ploton
    for jj = 1:num_Hard_Grasps
        figure(jj)
        subplot(4,1,1)
        plot(Data_Hard(jj).time, Data_Hard(jj).strain )
        xlabel('Index')
        ylabel('Strain ')

        subplot(4,1,2)
        plot(Data_Hard(jj).time, Data_Hard(jj).strainDot )
        xlabel('Index')
        ylabel('Straindot ')

        subplot(4,1,3)
        plot(Data_Hard(jj).time, Data_Hard(jj).strainDotDot )
        xlabel('Index')
        ylabel('Straindotdot ')

        subplot(4,1,4)
        plot(Data_Hard(jj).time, Data_Hard(jj).stress )
        xlabel('Index')
        ylabel('stress ')

        suptitle('Hard tissue')
    end


    for jj = 1:num_Soft_Grasps
        figure(jj+num_Hard_Grasps)
        subplot(4,1,1)
        plot(Data_Soft(jj).time, Data_Soft(jj).strain )
        xlabel('Index')
        ylabel('Strain ')

        subplot(4,1,2)
        plot(Data_Soft(jj).time, Data_Soft(jj).strainDot )
        xlabel('Index')
        ylabel('Straindot ')

        subplot(4,1,3)
        plot(Data_Soft(jj).time, Data_Soft(jj).strainDotDot )
        xlabel('Index')
        ylabel('Straindotdot ')

        subplot(4,1,4)
        plot(Data_Soft(jj).time, Data_Soft(jj).stress )
        xlabel('Index')
        ylabel('stress ')

        suptitle('Soft tissue')
    end

end



% data{1}.state = D_Hard.D1;
% data{2}.state = D_Soft.D1;
% data{1}.input = U_Hard.U1;
% data{2}.input = U_Soft.U1;
% 
% lambda = [0,0];
% 
% Parameters = DLS_TrainGeneral(data, lambda);

data_all = [];
data_all{1}.state = [];
data_all{1}.input = [];
data_all{2}.state = [];
data_all{2}.input = [];

%Get All Params for TLS
for ii = 1:num_Hard_Grasps
    
    data_all{1}.state = [data_all{1}.state; Data_Hard(ii).strainDotDot.*(Data_Hard(ii).firstTouch.^2), Data_Hard(ii).strainDot.*Data_Hard(ii).firstTouch, Data_Hard(ii).strain, Data_Hard(ii).strain.^2, Data_Hard(ii).strain.^3];
    data_all{1}.input = [data_all{1}.input; Data_Hard(ii).stress];

end
Params_Hard_All = pinv(data_all{1}.state)*data_all{1}.input;

x = lsqlin(data_all{1}.state,data_all{1}.input)
Residual_Hard=data_all{1}.input-data_all{1}.state*Params_Hard_All;


for ii = 1:num_Soft_Grasps
    data_all{2}.state = [data_all{2}.state; Data_Soft(ii).strainDotDot.*(Data_Soft(ii).firstTouch.^2), Data_Soft(ii).strainDot.*Data_Soft(ii).firstTouch, Data_Soft(ii).strain, Data_Soft(ii).strain.^2, Data_Soft(ii).strain.^3];
    data_all{2}.input  = [data_all{2}.input; Data_Soft(ii).stress];
end


Params_Soft_All = pinv(data_all{2}.state)*data_all{2}.input;

cov_soft = inv(data_all{2}.state'*data_all{2}.state);
%%
%DLS params
lambdas = 0:0.1:1;
for ind = 1:length(lambdas)
    %check dls params
    params_dls = DLS_TrainGeneral(data_all,[lambdas(ind),lambdas(ind)],0);
    
    Params_Hard_DLS(:,ind) = params_dls(:,1);
    Params_Soft_DLS(:,ind) = params_dls(:,2);
end

[nnh,orderh] = size(data_all{1}.state);
[nns,orders] = size(data_all{2}.state);

figure(1)
for ii = 1:orderh
    subplot(orderh,1,ii)
    plot(lambdas,Params_Hard_DLS(ii,:));
    xlabel('lambda');
    ylabel('parameter');
    str = sprintf('phi %d',ii);
    title(str);
end
suptitle('Hard tissue parameters')

figure(2)
for ii = 1:orders
    subplot(orders,1,ii)
    plot(lambdas,Params_Soft_DLS(ii,:));
    xlabel('lambda');
    ylabel('parameter');
    str = sprintf('phi %d',ii);
    title(str);
end
suptitle('Soft tissue parameters')


% Leave One out for Hard
for ii = 1:num_Hard_Grasps
    State_Hard_Train = [];
    Input_Hard_Train = [];
    State_Hard_Online = [];
    Input_Hard_Online = [];
    for jj = 1:num_Hard_Grasps
        if(ii ~= jj)
            State_Hard_Train = [State_Hard_Train; Data_Hard(jj).strainDotDot.*(Data_Hard(jj).firstTouch.^2), Data_Hard(jj).strainDot.*Data_Hard(jj).firstTouch, Data_Hard(jj).strain, Data_Hard(jj).strain.^2, Data_Hard(jj).strain.^3];
            Input_Hard_Train = [Input_Hard_Train; Data_Hard(jj).stress];
        else
            State_Hard_Online = [Data_Hard(jj).strainDotDot.*(Data_Hard(jj).firstTouch.^2), Data_Hard(jj).strainDot.*Data_Hard(jj).firstTouch, Data_Hard(jj).strain, Data_Hard(jj).strain.^2, Data_Hard(jj).strain.^3];
            Input_Hard_Online = [Data_Hard(jj).stress];
        end
    end
    %Classify?
    Params_Hard_Train(:,ii) = pinv(State_Hard_Train)*Input_Hard_Train;
    
    Error1 = sum(abs(Input_Hard_Online - State_Hard_Online*Params_Hard_Train(:,ii) )); %Hard
    Error2 = sum(abs(Input_Hard_Online - State_Hard_Online*Params_Soft_All )); %Soft
    if(Error1 < Error2)
        Classification_Hard(ii,:) = [Error1, Error2, hard_id, hard_id]; %Hard
    else
        Classification_Hard(ii,:) = [Error1, Error2, hard_id, soft_id]; %Soft
    end
end


% Leave One out for Soft
for ii = 1:num_Soft_Grasps
    State_Soft_Train = [];
    Input_Soft_Train = [];
    State_Soft_Online = [];
    Input_Soft_Online = [];
    for jj = 1:num_Soft_Grasps
        if(ii ~= jj)
            State_Soft_Train = [State_Soft_Train; Data_Soft(jj).strainDotDot.*(Data_Soft(jj).firstTouch.^2), Data_Soft(jj).strainDot, Data_Soft(jj).strain, Data_Soft(jj).strain.^2, Data_Soft(jj).strain.^3];
            Input_Soft_Train = [Input_Soft_Train; Data_Soft(jj).stress];
        else
            State_Soft_Online = [Data_Soft(jj).strainDotDot.*(Data_Soft(jj).firstTouch.^2), Data_Soft(jj).strainDot, Data_Soft(jj).strain, Data_Soft(jj).strain.^2, Data_Soft(jj).strain.^3];
            Input_Soft_Online = [Data_Soft(jj).stress];
        end
    end
    %Classify?
    Params_Soft_Train(:,ii) = pinv(State_Soft_Train)*Input_Soft_Train;
    
    Error1 = sum(abs(Input_Soft_Online - State_Soft_Online*Params_Hard_All )); %Hard
    Error2 = sum(abs(Input_Soft_Online - State_Soft_Online*Params_Soft_Train(:,ii) )); %Soft
    if(Error1 < Error2)
        Classification_Soft(ii,:) = [Error1, Error2, soft_id, hard_id]; %Hard
    else
        Classification_Soft(ii,:) = [Error1, Error2, soft_id, soft_id]; %Soft
    end
end









