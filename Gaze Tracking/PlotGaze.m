% import gaze data and plot
% Rod Dockter

clear all

% Go to M Drive (SLOOOOOOOOW)
oldFolder = cd('M:/Projects/GazeTracking/data/GazeData_12_8_14')

% list file names, must be matching edge for each gaze file
filenames_novice_EDGE = {'EDGE4_JJO_PegTransfer_Normal_12.08.2014-13.30.02.txt'; 'EDGE4_JJO_PegTransfer_TargetGaze_12.08.2014-13.37.27.txt'; 'EDGE4_JJO_PegTransfer_ToolGaze_12.08.2014-13.42.51.txt'  };
filenames_novice_GAZE = {'EDGE4_JJO_PegTransfer_Normal_GazeTrack.txt'; 'EDGE4_JJO_PegTransfer_TargetGaze_GazeTrack.txt'; 'EDGE4_JJO_PegTransfer_ToolGaze_GazeTrack.txt' };
filenames_novice_gazetype = {'Normal'; 'Target'; 'Tool'};
num_novices = length(filenames_novice_EDGE);

offsets_novice = [23,23,9]; % Amount Gaze video is ahead of edge, manually adjust via video

filenames_expert_EDGE = {'EDGE4_LP_PegTransfer_Normal2_12.08.2014-17.28.49.txt';  'EDGE4_LP_PegTransfer_TargetGaze_12.08.2014-17.33.25.txt'; 'EDGE4_LP_PegTransfer_ToolGaze_12.08.2014-17.36.28.txt'};
filenames_expert_GAZE = {'EDGE4_LP_PegTransfer_Normal2_GazeTrack.txt';  'EDGE4_LP_PegTransfer_TargetGaze_GazeTrack.txt'; 'EDGE4_LP_PegTransfer_ToolGaze_GazeTrack.txt'};
filenames_expert_gazetype = {'Normal'; 'Target'; 'Tool'};
num_experts = length(filenames_expert_EDGE);

offsets_expert = [7,8,6]; % Amount Gaze video is ahead of edge, manually adjust via video


%Load Data
for ii = 1 : num_novices
   DataNovice{ii}.edge = load(char(filenames_novice_EDGE(ii))); 
   DataNovice{ii}.gaze = load(char(filenames_novice_GAZE(ii)));
   DataNovice{ii}.type = filenames_novice_gazetype(ii);
   DataNovice{ii}.offset_time = offsets_novice(ii);
end

for ii = 1 : num_experts
   DataExpert{ii}.edge = load(char(filenames_expert_EDGE(ii))); 
   DataExpert{ii}.gaze = load(char(filenames_expert_GAZE(ii))); 
   DataExpert{ii}.type = filenames_expert_gazetype(ii);
   DataExpert{ii}.offset_time = offsets_expert(ii);
end

cd(oldFolder) 

%%

%Extract Data (time, J1, J2, gaze (x,y), and grasper force)
for ii = 1 : num_novices
   DataNovice{ii}.TimeEdge = DataNovice{ii}.edge(:,1);
   DataNovice{ii}.deltaTedge = mean(diff(DataNovice{ii}.TimeEdge));
   DataNovice{ii}.L.J1 = DataNovice{ii}.edge(:,2);
   DataNovice{ii}.L.J2 = DataNovice{ii}.edge(:,3);
   DataNovice{ii}.R.J1 = DataNovice{ii}.edge(:,8);
   DataNovice{ii}.R.J2 = DataNovice{ii}.edge(:,9);
   
   DataNovice{ii}.L.GrasperForce = DataNovice{ii}.edge(:,7);
   DataNovice{ii}.R.GrasperForce = DataNovice{ii}.edge(:,13);
   
   DataNovice{ii}.TimeGaze = ( DataNovice{ii}.gaze(:,2) - DataNovice{ii}.gaze(1,2) ) ./1000.0; %ms 2 s
   DataNovice{ii}.deltaTgaze = mean(diff(DataNovice{ii}.TimeGaze)); 
   DataNovice{ii}.X = DataNovice{ii}.gaze(:,3);
   DataNovice{ii}.Y = DataNovice{ii}.gaze(:,4);
end

for ii = 1 : num_experts
   DataExpert{ii}.TimeEdge = DataExpert{ii}.edge(:,1);
   DataExpert{ii}.deltaTedge = mean(diff(DataExpert{ii}.TimeEdge));
   DataExpert{ii}.L.J1 = DataExpert{ii}.edge(:,2);
   DataExpert{ii}.L.J2 = DataExpert{ii}.edge(:,3);
   DataExpert{ii}.R.J1 = DataExpert{ii}.edge(:,8);
   DataExpert{ii}.R.J2 = DataExpert{ii}.edge(:,9);
   
   DataExpert{ii}.L.GrasperForce = DataExpert{ii}.edge(:,7);
   DataExpert{ii}.R.GrasperForce = DataExpert{ii}.edge(:,13);
   
   DataExpert{ii}.TimeGaze = ( DataExpert{ii}.gaze(:,2)  - DataExpert{ii}.gaze(1,2) ) ./1000.0; %ms 2 s
   DataExpert{ii}.deltaTgaze = mean(diff(DataExpert{ii}.TimeGaze));
   DataExpert{ii}.X = DataExpert{ii}.gaze(:,3);
   DataExpert{ii}.Y = DataExpert{ii}.gaze(:,4);
end

%%

% Calculate Velocity
for ii = 1 : num_novices
    %tool
    DataNovice{ii}.L.J1dot = Calculate_velocity( DataNovice{ii}.L.J1, DataNovice{ii}.deltaTedge, 'holobrodko'); 
    DataNovice{ii}.L.J2dot = Calculate_velocity( DataNovice{ii}.L.J2, DataNovice{ii}.deltaTedge, 'holobrodko');
    DataNovice{ii}.R.J1dot = Calculate_velocity( DataNovice{ii}.R.J1, DataNovice{ii}.deltaTedge, 'holobrodko');
    DataNovice{ii}.R.J2dot = Calculate_velocity( DataNovice{ii}.R.J2, DataNovice{ii}.deltaTedge, 'holobrodko');
    DataNovice{ii}.L.ToolVel = sqrt( DataNovice{ii}.L.J1dot.^2 + DataNovice{ii}.L.J2dot.^2 );
    DataNovice{ii}.R.ToolVel = sqrt( DataNovice{ii}.R.J1dot.^2 + DataNovice{ii}.R.J2dot.^2 );
    DataNovice{ii}.L.ToolVel_scaled = DataNovice{ii}.L.ToolVel ./ max(DataNovice{ii}.L.ToolVel);
    DataNovice{ii}.R.ToolVel_scaled = DataNovice{ii}.R.ToolVel ./ max(DataNovice{ii}.R.ToolVel);
    
    %gaze
    DataNovice{ii}.Xdot = Calculate_velocity( DataNovice{ii}.X, DataNovice{ii}.deltaTgaze, 'holobrodko');
    DataNovice{ii}.Ydot = Calculate_velocity( DataNovice{ii}.Y, DataNovice{ii}.deltaTgaze, 'holobrodko');
    DataNovice{ii}.GazeVel = sqrt( DataNovice{ii}.Xdot.^2 + DataNovice{ii}.Ydot.^2 );
    DataNovice{ii}.GazeVel_scaled = DataNovice{ii}.GazeVel ./ max(DataNovice{ii}.GazeVel);
end

for ii = 1 : num_experts
    %tool
    DataExpert{ii}.L.J1dot = Calculate_velocity( DataExpert{ii}.L.J1, DataExpert{ii}.deltaTedge, 'holobrodko');
    DataExpert{ii}.L.J2dot = Calculate_velocity( DataExpert{ii}.L.J2, DataExpert{ii}.deltaTedge, 'holobrodko'); 
    DataExpert{ii}.R.J1dot = Calculate_velocity( DataExpert{ii}.R.J1, DataExpert{ii}.deltaTedge, 'holobrodko'); 
    DataExpert{ii}.R.J2dot = Calculate_velocity( DataExpert{ii}.R.J2, DataExpert{ii}.deltaTedge, 'holobrodko');
    DataExpert{ii}.L.ToolVel = sqrt( DataExpert{ii}.L.J1dot.^2 + DataExpert{ii}.L.J2dot.^2 );
    DataExpert{ii}.R.ToolVel = sqrt( DataExpert{ii}.R.J1dot.^2 + DataExpert{ii}.R.J2dot.^2 );
    DataExpert{ii}.L.ToolVel_scaled = DataExpert{ii}.L.ToolVel ./ max(DataExpert{ii}.L.ToolVel);
    DataExpert{ii}.R.ToolVel_scaled = DataExpert{ii}.R.ToolVel ./ max(DataExpert{ii}.R.ToolVel);
    
    %gaze
    DataExpert{ii}.Xdot = Calculate_velocity( DataExpert{ii}.X, DataExpert{ii}.deltaTgaze, 'holobrodko');
    DataExpert{ii}.Ydot = Calculate_velocity( DataExpert{ii}.Y, DataExpert{ii}.deltaTgaze, 'holobrodko');
    DataExpert{ii}.GazeVel = sqrt( DataExpert{ii}.Xdot.^2 + DataExpert{ii}.Ydot.^2 );
    DataExpert{ii}.GazeVel_scaled = DataExpert{ii}.GazeVel ./ max(DataExpert{ii}.GazeVel);
end

%%

%geting plot indices

start_time = 2; %start after 1 sec
end_time = 45;

for ii = 1 : num_novices
    %indexes
    DataNovice{ii}.index_start_edge = find(DataNovice{ii}.TimeEdge > start_time, 1, 'first'); 
    DataNovice{ii}.index_end_edge = find(DataNovice{ii}.TimeEdge > end_time, 1, 'first'); 
    DataNovice{ii}.index_start_gaze = find(DataNovice{ii}.TimeGaze > start_time + DataNovice{ii}.offset_time, 1, 'first'); 
    DataNovice{ii}.index_end_gaze = find(DataNovice{ii}.TimeGaze > end_time + DataNovice{ii}.offset_time, 1, 'first'); 
    
end

for ii = 1 : num_experts
    %indexes
    DataExpert{ii}.index_start_edge = find(DataExpert{ii}.TimeEdge > start_time, 1, 'first'); 
    DataExpert{ii}.index_end_edge = find(DataExpert{ii}.TimeEdge > end_time, 1, 'first'); 
    DataExpert{ii}.index_start_gaze = find(DataExpert{ii}.TimeGaze > start_time + DataExpert{ii}.offset_time, 1, 'first'); 
    DataExpert{ii}.index_end_gaze = find(DataExpert{ii}.TimeGaze > end_time + DataExpert{ii}.offset_time, 1, 'first'); 
end


%%

%Plotting
for ii = 1 : num_novices
    figure(ii)
    subplot(2,1,1);
    plot(DataNovice{ii}.TimeGaze(DataNovice{ii}.index_start_gaze:DataNovice{ii}.index_end_gaze) - DataNovice{ii}.offset_time,DataNovice{ii}.GazeVel_scaled(DataNovice{ii}.index_start_gaze:DataNovice{ii}.index_end_gaze),'-','MarkerSize',2,'MarkerFaceColor',[1 0 0],'Color',[1 0 0])
    hold on
    plot(DataNovice{ii}.TimeEdge(DataNovice{ii}.index_start_edge:DataNovice{ii}.index_end_edge),DataNovice{ii}.L.ToolVel_scaled(DataNovice{ii}.index_start_edge:DataNovice{ii}.index_end_edge),'+','MarkerSize',2,'MarkerFaceColor',[0 1 0],'Color',[0 1 0])
    hold on
    plot(DataNovice{ii}.TimeEdge(DataNovice{ii}.index_start_edge:DataNovice{ii}.index_end_edge),DataNovice{ii}.R.ToolVel_scaled(DataNovice{ii}.index_start_edge:DataNovice{ii}.index_end_edge),'x','MarkerSize',2,'MarkerFaceColor',[0 0 1],'Color',[0 0 1])
    hold off
    
    titler = strcat('Novice Gaze vs. Tool. Gaze Type = ',DataNovice{ii}.type);
    title(titler)
    xlabel('Time (Seconds)')
    ylabel('Velocity (Scaled)')
    legend('Gaze Velocity','Left Tool Velocity','Right Tool Velocity','Location','BestOutside')

    subplot(2,1,2);
    plot(DataNovice{ii}.TimeEdge(DataNovice{ii}.index_start_edge:DataNovice{ii}.index_end_edge),DataNovice{ii}.L.GrasperForce(DataNovice{ii}.index_start_edge:DataNovice{ii}.index_end_edge),'o','MarkerSize',2,'MarkerFaceColor',[1 0 0],'Color',[1 0 0])
    hold on
    plot(DataNovice{ii}.TimeEdge(DataNovice{ii}.index_start_edge:DataNovice{ii}.index_end_edge),DataNovice{ii}.R.GrasperForce(DataNovice{ii}.index_start_edge:DataNovice{ii}.index_end_edge),'-','MarkerSize',2,'MarkerFaceColor',[0 1 0],'Color',[0 1 0])
    hold off
    
    title('Grasping Force')
    xlabel('Time (Seconds)')
    ylabel('Grasper Force (Scaled)')
    legend('Left Grasper','Right Grasper','Location','BestOutside')

    
end


for ii = 1 : num_experts
    figure(ii + num_novices)
    subplot(2,1,1);
    plot(DataExpert{ii}.TimeGaze(DataExpert{ii}.index_start_gaze:DataExpert{ii}.index_end_gaze) - DataExpert{ii}.offset_time,DataExpert{ii}.GazeVel_scaled(DataExpert{ii}.index_start_gaze:DataExpert{ii}.index_end_gaze),'-','MarkerSize',2,'MarkerFaceColor',[1 0 0],'Color',[1 0 0])
    hold on
    plot(DataExpert{ii}.TimeEdge(DataExpert{ii}.index_start_edge:DataExpert{ii}.index_end_edge),DataExpert{ii}.L.ToolVel_scaled(DataExpert{ii}.index_start_edge:DataExpert{ii}.index_end_edge),'+','MarkerSize',2,'MarkerFaceColor',[0 1 0],'Color',[0 1 0])
    hold on
    plot(DataExpert{ii}.TimeEdge(DataExpert{ii}.index_start_edge:DataExpert{ii}.index_end_edge),DataExpert{ii}.R.ToolVel_scaled(DataExpert{ii}.index_start_edge:DataExpert{ii}.index_end_edge),'x','MarkerSize',2,'MarkerFaceColor',[0 0 1],'Color',[0 0 1])
    hold off
    
    titler = strcat('Expert Gaze vs. Tool. Gaze Type = ',DataExpert{ii}.type);
    title(titler)
    xlabel('Time (Seconds)')
    ylabel('Velocity (Scaled)')
    legend('Gaze Velocity','Left Tool Velocity','Right Tool Velocity','Location','BestOutside')

    subplot(2,1,2);
    plot(DataExpert{ii}.TimeEdge(DataExpert{ii}.index_start_edge:DataExpert{ii}.index_end_edge),DataExpert{ii}.L.GrasperForce(DataExpert{ii}.index_start_edge:DataExpert{ii}.index_end_edge),'o','MarkerSize',2,'MarkerFaceColor',[1 0 0],'Color',[1 0 0])
    hold on
    plot(DataExpert{ii}.TimeEdge(DataExpert{ii}.index_start_edge:DataExpert{ii}.index_end_edge),DataExpert{ii}.R.GrasperForce(DataExpert{ii}.index_start_edge:DataExpert{ii}.index_end_edge),'-','MarkerSize',2,'MarkerFaceColor',[0 1 0],'Color',[0 1 0])
    hold off
    
    title('Grasping Force')
    xlabel('Time (Seconds)')
    ylabel('Grasper Force (Scaled)')
    legend('Left Grasper','Right Grasper','Location','BestOutside')

    
end


















