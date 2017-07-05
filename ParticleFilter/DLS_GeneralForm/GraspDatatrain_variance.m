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


% Do diffs
for jj = 1:num_Hard_Grasps

    dt = 0.001;
    
    %The default span for the moving average is 5.
    Data_Hard(jj).strainDot = smooth(diff(Data_Hard(jj).strain)/dt ,7);
    Data_Hard(jj).strainDot = [Data_Hard(jj).strainDot(1); Data_Hard(jj).strainDot];

    Data_Hard(jj).strainDotDot = smooth(diff(Data_Hard(jj).strainDot)/dt,7);
    Data_Hard(jj).strainDotDot = [Data_Hard(jj).strainDotDot(1); Data_Hard(jj).strainDotDot];
end


for jj = 1:num_Soft_Grasps
    dt = 0.001;

    %The default span for the moving average is 5.
    Data_Soft(jj).strainDot = smooth(diff(Data_Soft(jj).strain)/dt,7);
    Data_Soft(jj).strainDot = [Data_Soft(jj).strainDot(1); Data_Soft(jj).strainDot];

    Data_Soft(jj).strainDotDot = smooth(diff(Data_Soft(jj).strainDot)/dt,7);
    Data_Soft(jj).strainDotDot = [Data_Soft(jj).strainDotDot(1); Data_Soft(jj).strainDotDot];
end

%Get All Params for TLS
for ii = 1:num_Hard_Grasps
    
    Data_temp = [Data_Hard(ii).strainDotDot.*(Data_Hard(ii).firstTouch.^2), Data_Hard(ii).strainDot.*Data_Hard(ii).firstTouch, ones(length(Data_Hard(ii).strain),1), Data_Hard(ii).strain, Data_Hard(ii).strain.^2, Data_Hard(ii).strain.^3];
    Input_temp = Data_Hard(ii).stress / 1000;
    
    Params_Hard_LOO(:,ii) = pinv(Data_temp)*Input_temp;
    
    data_all{1}.state = [data_all{1}.state; Data_temp];
    data_all{1}.input = [data_all{1}.input; Input_temp];

end

Params_Hard_All = pinv(data_all{1}.state)*data_all{1}.input
Params_Hard_std = std(Params_Hard_LOO')


X = data_all{1}.state(:,3);
Y = data_all{1}.state(:,4);
Z = data_all{1}.state(:,5);


figure
j=1;
subplot(5,1,j)
plot(Params_Hard_LOO(j,:));
j=j+1;
subplot(5,1,j)
plot(Params_Hard_LOO(j,:));
j=j+1;
subplot(5,1,j)
plot(Params_Hard_LOO(j,:));
j=j+1;
subplot(5,1,j)
plot(Params_Hard_LOO(j,:));
j=j+1;
subplot(5,1,j)
plot(Params_Hard_LOO(j,:));
suptitle('parameters variation over grasps (hard)')



for ii = 1:num_Soft_Grasps
    
    
    Data_temp = [Data_Soft(ii).strainDotDot.*(Data_Soft(ii).firstTouch.^2), Data_Soft(ii).strainDot.*Data_Soft(ii).firstTouch, ones(length(Data_Soft(ii).strain),1), Data_Soft(ii).strain, Data_Soft(ii).strain.^2, Data_Soft(ii).strain.^3];
    Input_temp = Data_Soft(ii).stress / 1000;
    
    Params_Soft_LOO(:,ii) = pinv(Data_temp)*Input_temp;
    
    data_all{2}.state = [data_all{2}.state; Data_temp];
    data_all{2}.input  = [data_all{2}.input; Input_temp];
end


Params_Soft_All = pinv(data_all{2}.state)*data_all{2}.input
Params_Soft_std = std(Params_Soft_LOO')

figure
j=1;
subplot(5,1,j)
plot(Params_Soft_LOO(j,:));
j=j+1;
subplot(5,1,j)
plot(Params_Soft_LOO(j,:));
j=j+1;
subplot(5,1,j)
plot(Params_Soft_LOO(j,:));
j=j+1;
subplot(5,1,j)
plot(Params_Soft_LOO(j,:));
j=j+1;
subplot(5,1,j)
plot(Params_Soft_LOO(j,:));
suptitle('parameters variation over grasps (soft)')
