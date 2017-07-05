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


data_all = [];
data_all{1}.state = [];
data_all{1}.input = [];
data_all{2}.state = [];
data_all{2}.input = [];

%Get All Params for TLS
for ii = 1:num_Hard_Grasps
    
    data_all{1}.state = [data_all{1}.state; Data_Hard(ii).strain, Data_Hard(ii).strain.^2, Data_Hard(ii).strain.^3];
    data_all{1}.input = [data_all{1}.input; Data_Hard(ii).stress];

end
Params_Hard_All = pinv(data_all{1}.state)*data_all{1}.input;

for ii = 1:num_Soft_Grasps
    data_all{2}.state = [data_all{2}.state; Data_Soft(ii).strain, Data_Soft(ii).strain.^2, Data_Soft(ii).strain.^3];
    data_all{2}.input  = [data_all{2}.input; Data_Soft(ii).stress];
end
Params_Soft_All = pinv(data_all{2}.state)*data_all{2}.input;


%DLS params
lambdas = 0:0.1:1;
for ind = 1:length(lambdas)
    %check dls params
    params_dls = DLS_TrainGeneral(data_all,[lambdas(ind),lambdas(ind)],0);
    
    Params_Hard_DLS(:,ind) = params_dls(:,1);
    Params_Soft_DLS(:,ind) = params_dls(:,2);
end

%Get order of each data
[nnh,orderh] = size(data_all{1}.state);
[nns,orders] = size(data_all{2}.state);


%Plot change in parameters as lambda increases
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


data_train = [];
data_train{1}.state = data_all{1}.state;
data_train{1}.input = data_all{1}.input;
data_train{2}.state = data_all{2}.state;
data_train{2}.input = data_all{2}.input;

% Leave One out for Hard
for ii = 1:num_Hard_Grasps
    data_train{1}.state = [];
    data_train{1}.input = [];
    State_Hard_Online = [];
    Input_Hard_Online = [];
    for jj = 1:num_Hard_Grasps
        if(ii ~= jj)
            data_train{1}.state = [data_train{1}.state; Data_Hard(jj).strain, Data_Hard(jj).strain.^2, Data_Hard(jj).strain.^3];
            data_train{1}.input = [data_train{1}.input; Data_Hard(jj).stress];
        else
            State_Hard_Online = [Data_Hard(jj).strain, Data_Hard(jj).strain.^2, Data_Hard(jj).strain.^3];
            Input_Hard_Online = [Data_Hard(jj).stress];
        end
    end
    %TLS parameters
    Params_Hard_Train(:,ii) = pinv(data_train{1}.state)*data_train{1}.input;
    
    %TLS classify
    Error1 = sum(abs(Input_Hard_Online - State_Hard_Online*Params_Hard_Train(:,ii) )); %Hard
    Error2 = sum(abs(Input_Hard_Online - State_Hard_Online*Params_Soft_All )); %Soft
    if(Error1 < Error2)
        Classification_Hard(ii,:) = [Error1, Error2, hard_id, hard_id]; %Hard
    else
        Classification_Hard(ii,:) = [Error1, Error2, hard_id, soft_id]; %Soft
    end
    
    
    %DLS parameters
    for ind = 1:length(lambdas)
        params_dls = DLS_TrainGeneral(data_train,[lambdas(ind),lambdas(ind)],0);
        Params_Hard_Train_DLS(:,ind) = params_dls(:,1);
        
        Error1 = sum(abs(Input_Hard_Online - State_Hard_Online*Params_Hard_Train_DLS(:,ind) )); %Hard
        Error2 = sum(abs(Input_Hard_Online - State_Hard_Online*Params_Soft_DLS(:,ind) )); %Soft
        if(Error1 < Error2)
            Classification_DLS_Hard{ind}(ii,:) = [Error1, Error2, hard_id, hard_id]; %Hard
        else
            Classification_DLS_Hard{ind}(ii,:) = [Error1, Error2, hard_id, soft_id]; %Soft
        end
    end
    

end


% Leave One out for Soft
for ii = 1:num_Soft_Grasps
    data_train{2}.state = [];
    data_train{2}.input = [];
    State_Soft_Online = [];
    Input_Soft_Online = [];
    for jj = 1:num_Soft_Grasps
        if(ii ~= jj)
            data_train{2}.state = [data_train{2}.state; Data_Soft(jj).strain, Data_Soft(jj).strain.^2, Data_Soft(jj).strain.^3];
            data_train{2}.input = [data_train{2}.input; Data_Soft(jj).stress];
        else
            State_Soft_Online = [Data_Soft(jj).strain, Data_Soft(jj).strain.^2, Data_Soft(jj).strain.^3];
            Input_Soft_Online = [Data_Soft(jj).stress];
        end
    end
    %TLS parameters
    Params_Soft_Train(:,ii) = pinv(data_train{2}.state)*data_train{2}.input;
    
    %TLS classify
    Error1 = sum(abs(Input_Soft_Online - State_Soft_Online*Params_Hard_All )); %Hard
    Error2 = sum(abs(Input_Soft_Online - State_Soft_Online*Params_Soft_Train(:,ii) )); %Soft
    if(Error1 < Error2)
        Classification_Soft(ii,:) = [Error1, Error2, soft_id, hard_id]; %Hard
    else
        Classification_Soft(ii,:) = [Error1, Error2, soft_id, soft_id]; %Soft
    end
    
    
    %DLS parameters
    for ind = 1:length(lambdas)
        params_dls = DLS_TrainGeneral(data_train,[lambdas(ind),lambdas(ind)],0);
        Params_Soft_Train_DLS(:,ind) = params_dls(:,2);
        
        Error1 = sum(abs(Input_Soft_Online - State_Soft_Online*Params_Hard_DLS(:,ind) )); %Hard
        Error2 = sum(abs(Input_Soft_Online - State_Soft_Online*Params_Soft_Train_DLS(:,ind) )); %Soft
        if(Error1 < Error2)
            Classification_DLS_Soft{ind}(ii,:) = [Error1, Error2, soft_id, hard_id]; %Hard
        else
            Classification_DLS_Soft{ind}(ii,:) = [Error1, Error2, soft_id, soft_id]; %Soft
        end
    end
end









