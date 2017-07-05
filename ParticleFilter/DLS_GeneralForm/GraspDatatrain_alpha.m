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
    data_all{1}.state = [data_all{1}.state; Data_Hard(ii).strainDotDot.*(Data_Hard(ii).firstTouch.^2), Data_Hard(ii).strainDot.*Data_Hard(ii).firstTouch, ones(length(Data_Hard(ii).strain),1), Data_Hard(ii).strain, Data_Hard(ii).strain.^2, Data_Hard(ii).strain.^3];
    data_all{1}.input = [data_all{1}.input; Data_Hard(ii).stress / 1000];

end
Params_Hard_All = pinv(data_all{1}.state)*data_all{1}.input;

for ii = 1:num_Soft_Grasps
    data_all{2}.state = [data_all{2}.state; Data_Soft(ii).strainDotDot.*(Data_Soft(ii).firstTouch.^2), Data_Soft(ii).strainDot.*Data_Soft(ii).firstTouch, ones(length(Data_Soft(ii).strain),1), Data_Soft(ii).strain, Data_Soft(ii).strain.^2, Data_Soft(ii).strain.^3];
    data_all{2}.input  = [data_all{2}.input; Data_Soft(ii).stress / 1000];
end
Params_Soft_All = pinv(data_all{2}.state)*data_all{2}.input;



% Leave One out for Hard
for ii = 1:num_Hard_Grasps
    State_Hard_Train = [];
    Input_Hard_Train = [];
    State_Hard_Online = [];
    Input_Hard_Online = [];
    for jj = 1:num_Hard_Grasps
        if(ii ~= jj)
            State_Hard_Train = [State_Hard_Train; Data_Hard(jj).strainDotDot.*(Data_Hard(jj).firstTouch.^2), Data_Hard(jj).strainDot.*Data_Hard(jj).firstTouch, ones(length(Data_Hard(jj).strain),1), Data_Hard(jj).strain, Data_Hard(jj).strain.^2, Data_Hard(jj).strain.^3];
            Input_Hard_Train = [Input_Hard_Train; Data_Hard(jj).stress / 1000];
        else
            State_Hard_Online = [Data_Hard(jj).strainDotDot.*(Data_Hard(jj).firstTouch.^2), Data_Hard(jj).strainDot.*Data_Hard(jj).firstTouch, ones(length(Data_Hard(jj).strain),1), Data_Hard(jj).strain, Data_Hard(jj).strain.^2, Data_Hard(jj).strain.^3];
            Input_Hard_Online = [Data_Hard(jj).stress / 1000];
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
    
    Error_CumSum_Hard{ii} = [cumsum(abs(Input_Hard_Online - State_Hard_Online*Params_Hard_Train(:,ii) )), cumsum(abs(Input_Hard_Online - State_Hard_Online*Params_Soft_All ))] ;
end


% Leave One out for Soft
for ii = 1:num_Soft_Grasps
    State_Soft_Train = [];
    Input_Soft_Train = [];
    State_Soft_Online = [];
    Input_Soft_Online = [];
    for jj = 1:num_Soft_Grasps
        if(ii ~= jj)
            State_Soft_Train = [State_Soft_Train; Data_Soft(jj).strainDotDot.*(Data_Soft(jj).firstTouch.^2), Data_Soft(jj).strainDot.*Data_Soft(jj).firstTouch, ones(length(Data_Soft(jj).strain),1), Data_Soft(jj).strain, Data_Soft(jj).strain.^2, Data_Soft(jj).strain.^3];
            Input_Soft_Train = [Input_Soft_Train; Data_Soft(jj).stress / 1000];
        else
            State_Soft_Online = [Data_Soft(jj).strainDotDot.*(Data_Soft(jj).firstTouch.^2), Data_Soft(jj).strainDot.*Data_Soft(jj).firstTouch, ones(length(Data_Soft(jj).strain),1), Data_Soft(jj).strain, Data_Soft(jj).strain.^2, Data_Soft(jj).strain.^3];
            Input_Soft_Online = [Data_Soft(jj).stress / 1000];
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
    
    Error_CumSum_Soft{ii} = [cumsum(abs(Input_Soft_Online - State_Soft_Online*Params_Hard_All )), cumsum(abs(Input_Soft_Online - State_Soft_Online*Params_Soft_Train(:,ii) ))] ;
end


%Plot error sums to figure out alpha
figure
for ii = 1:num_Hard_Grasps
    h1 = plot(1:length(Error_CumSum_Hard{ii}(:,1)),Error_CumSum_Hard{ii}(:,1),'r-');
    h2 = plot(1:length(Error_CumSum_Hard{ii}(:,2)),Error_CumSum_Hard{ii}(:,2),'b-');
    hold on
end
hold off
xlabel('index')
ylabel('error')
title('change in error over time (Hard)')
legend([h1,h2],'Hard','Soft')


figure
for ii = 1:num_Soft_Grasps
    h1 = plot(1:length(Error_CumSum_Soft{ii}(:,1)),Error_CumSum_Soft{ii}(:,1),'r-');
    h2 = plot(1:length(Error_CumSum_Soft{ii}(:,2)),Error_CumSum_Soft{ii}(:,2),'b-');
    hold on
end
hold off
xlabel('index')
ylabel('error')
title('change in error over time (Soft)')
legend([h1,h2],'Hard','Soft')


alpha_thresh = 0.2;

ALPHA_SCALE_H = 100;
ALPHA_SCALE_S = 100;

FIG_W = 3.5;     % Width of actual figure  
FIG_H = 2;     % Height of actual figure
FIG_UNITS = 'inches'; % units for W&H
FIG_RES = 400; % figure resolution in dpi

%Plot alphas
figure
for ii = 1:num_Hard_Grasps
    diffe = (Error_CumSum_Hard{ii}(:,2) - Error_CumSum_Hard{ii}(:,1)) ./ (Error_CumSum_Hard{ii}(:,2)+ALPHA_SCALE_H);
    h1 = plot([1:length(diffe)]*0.001,diffe,'r-');
    idx_hard(ii) = find(diffe > alpha_thresh, 1, 'first');
    mean_error_hard(ii) = mean(Error_CumSum_Hard{ii}(:,1));
    mean_diff_hard(ii) = mean(Error_CumSum_Hard{ii}(:,2) - Error_CumSum_Hard{ii}(:,1));
    hold on
end

for ii = 1:num_Soft_Grasps
    diffe = (Error_CumSum_Soft{ii}(:,1) - Error_CumSum_Soft{ii}(:,2)) ./ (Error_CumSum_Soft{ii}(:,1)+ALPHA_SCALE_S);
    h2 = plot([1:length(diffe)]*0.001,diffe,'b-');
    idx_soft(ii) = find(diffe > alpha_thresh, 1, 'first');
    mean_error_soft(ii) = mean(Error_CumSum_Soft{ii}(:,2));
    mean_diff_soft(ii) = mean(Error_CumSum_Soft{ii}(:,1) - Error_CumSum_Soft{ii}(:,2));
    hold on
end
hold off
axis([0,0.11,0,1])
xlabel('Time (s)','FontSize',10,'FontName','Times')
ylabel('Alpha','FontSize',10,'FontName','Times')
legend([h1(1),h2(1)],{'Hard','Soft'},'FontSize',8,'Location','southeast','FontName','Times')

% set it's W and H w/o messing up the position on the screen
set(gcf,'PaperPositionMode','auto', 'units', FIG_UNITS)
FIG_SZ = get(gcf, 'position');
FIG_SZ(3:end) = [FIG_W FIG_H];
set(gcf, 'position', FIG_SZ);


% Save the figure to file
print(gcf, 'C:\Users\MRDLAB\Documents\Rod Dockter\Surgical Robotics\mrdPublications\ICRA_Letters_Grasper\Images\alpha_combined', '-dpng', ['-r' num2str(FIG_RES)]) % Raster



max(idx_hard)
max(idx_soft)

mean(mean_error_hard)
mean(mean_error_soft)

mean(mean_diff_hard)
mean(mean_diff_soft)
