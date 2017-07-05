d = load ('MSIM_2015_6_16_19_41_44.txt');

[mimic.steps,mimic.inputs] = size(d)

fignum = 1;
fontS = 14;

%ONE LOOK UP TABLE TO RULE THEM AAAAAALLLLLLLLL
%Lookup_Table = {}
%Lookup_Table[1] = "Pick and Place Jacks"
%Lookup_Table[2] = "Pick and Place - Bowls"
%Lookup_Table[3] = "Background Table/Board"
%Lookup_Table[4] = "The Other - Tools"
%Lookup_Table[5] = "Generic Task Rings"
%Lookup_Table[6] = "Peg Board - Pegs"
%Lookup_Table[7] = "Peg Board - Special Pegs"
%Lookup_Table[8] = "Match Board - Shapes (numbers, letters)"
%Lookup_Table[9] = "Match Board - Goals"
%Lookup_Table[10] = "Match Board - Doors"
%Lookup_Table[11] = "Match Board - Sliders"
%Lookup_Table[12] = "Ring and Rail - Rail"
%Lookup_Table[13] = "Ring and Rail - Target Disc"
%Lookup_Table[14] = "Rope Walk - Rope"  
%Lookup_Table[15] = "Rope Walk - Hook"
%Lookup_Table[16] = "Background Tissue"
%Lookup_Table[17] = "Generic Needle"
%Lookup_Table[18] = "Generic Suture Thread"
%Lookup_Table[19] = "Thread The Rings - Goals Rings"
%Lookup_Table[20] = "Thread The Rings - Ring Stands"
%Lookup_Table[21] = "Needle Targeting - Target Pads"
%Lookup_Table[22] = "Needle Targeting - Needle Rack"

%pedals: 1=energy1s, 2=swap, 3=energy2p, 4=lclutch, 5=energy1p, 6=camera,7=rclutch, 8=clutch, 9=energy2s

%%

%time
mimic.t = d(:,1);
mimic.deltaT = diff(mimic.t);

mimic.timeStep = mean(mimic.deltaT);

%%%%%%%%%%%%%%%%%positions slave
sp_index = 2;
sp_total = 6;
mimic.tool{1}.slave.X = d(:,sp_index);
mimic.tool{1}.slave.Y = d(:,sp_index+1);
mimic.tool{1}.slave.Z = d(:,sp_index+2);
mimic.tool{2}.slave.X = d(:,sp_index+3);
mimic.tool{2}.slave.Y = d(:,sp_index+4);
mimic.tool{2}.slave.Z = d(:,sp_index+5);

%%%%%%%%%%%%%%%%%%%%%%%%orientations (quaternion)
so_index = sp_index + sp_total;
so_total = 8;

mimic.tool{1}.slave.q0 = d(:,so_index);
mimic.tool{1}.slave.q1 = d(:,so_index+1);
mimic.tool{1}.slave.q2 = d(:,so_index+2);
mimic.tool{1}.slave.q3 = d(:,so_index+3);

mimic.tool{2}.slave.q0 = d(:,so_index+4);
mimic.tool{2}.slave.q1 = d(:,so_index+5);
mimic.tool{2}.slave.q2 = d(:,so_index+6);
mimic.tool{2}.slave.q3 = d(:,so_index+7);


%%%%%%%%%%%%%Master input data
m_index = so_index + so_total;
m_total = 17;

%positions master
mimic.tool{1}.master.X = d(:,m_index);
mimic.tool{1}.master.Y = d(:,m_index+1);
mimic.tool{1}.master.Z = d(:,m_index+2);
mimic.tool{2}.master.X = d(:,m_index+3);
mimic.tool{2}.master.Y = d(:,m_index+4);
mimic.tool{2}.master.Z = d(:,m_index+5);

%orientations (master)
mimic.tool{1}.master.q0 = d(:,m_index+6);
mimic.tool{1}.master.q1 = d(:,m_index+7);
mimic.tool{1}.master.q2 = d(:,m_index+8);
mimic.tool{1}.master.q3 = d(:,m_index+9);

mimic.tool{2}.master.q0 = d(:,m_index+10);
mimic.tool{2}.master.q1 = d(:,m_index+11);
mimic.tool{2}.master.q2 = d(:,m_index+12);
mimic.tool{2}.master.q3 = d(:,m_index+13);

mimic.tool{1}.master.grasp = d(:,m_index+14);
mimic.tool{2}.master.grasp = d(:,m_index+15);

mimic.pedals = d(:,m_index+16);


%%%%%%%%%%%%%%%%%%%%%%%%interaction touch data
to_index = m_index + m_total;
to_total = 12;

mimic.is_touch = d(:,to_index);

mimic.tool{1}.slave.jaw_interact = d(:,to_index+1);
mimic.tool{1}.slave.aux_interact = d(:,to_index+2);
mimic.tool{1}.slave.normal_force = d(:,to_index+3);
mimic.tool{1}.slave.shear_force = d(:,to_index+4);
mimic.tool{2}.slave.jaw_interact = d(:,to_index+5);
mimic.tool{2}.slave.aux_interact = d(:,to_index+6);
mimic.tool{2}.slave.normal_force = d(:,to_index+7);
mimic.tool{2}.slave.shear_force = d(:,to_index+8);
mimic.tool{3}.slave.jaw_interact = d(:,to_index+9);
mimic.tool{3}.slave.normal_force = d(:,to_index+10);
mimic.tool{3}.slave.shear_force = d(:,to_index+11);

%check on total columns
totalcols = to_index + to_total - 1


%get velocity and acceleration
mimic.tool{1}.slave.dX = Calculate_velocity( mimic.tool{1}.slave.X, mimic.timeStep, 'holobrodko');
mimic.tool{1}.slave.dY = Calculate_velocity( mimic.tool{1}.slave.Y, mimic.timeStep, 'holobrodko');
mimic.tool{1}.slave.dZ = Calculate_velocity( mimic.tool{1}.slave.Z, mimic.timeStep, 'holobrodko');
mimic.tool{2}.slave.dX = Calculate_velocity( mimic.tool{2}.slave.X, mimic.timeStep, 'holobrodko');
mimic.tool{2}.slave.dY = Calculate_velocity( mimic.tool{2}.slave.Y, mimic.timeStep, 'holobrodko');
mimic.tool{2}.slave.dZ = Calculate_velocity( mimic.tool{2}.slave.Z, mimic.timeStep, 'holobrodko');



%convert quaternions to euler angles for each tool
mimic.tool{1}.slave.EA = SpinCalc('QtoEA123',[mimic.tool{1}.slave.q0,mimic.tool{1}.slave.q1,mimic.tool{1}.slave.q2,mimic.tool{1}.slave.q3],0.05,1) .* (pi/180);
mimic.tool{2}.slave.EA = SpinCalc('QtoEA123',[mimic.tool{2}.slave.q0,mimic.tool{2}.slave.q1,mimic.tool{2}.slave.q2,mimic.tool{2}.slave.q3],0.05,1) .* (pi/180);

mimic.tool{1}.master.EA = SpinCalc('QtoEA123',[mimic.tool{1}.master.q0,mimic.tool{1}.master.q1,mimic.tool{1}.master.q2,mimic.tool{1}.master.q3],0.05,1) .* (pi/180);
mimic.tool{2}.master.EA = SpinCalc('QtoEA123',[mimic.tool{2}.master.q0,mimic.tool{2}.master.q1,mimic.tool{2}.master.q2,mimic.tool{2}.master.q3],0.05,1) .* (pi/180);


%% grasping (0.6 = open, -0.6 = closed) probably radians

fignum = fignum + 1;
figure(fignum)
scatter(mimic.t,mimic.tool{1}.master.grasp)
xlabel('time (s)')
ylabel('grap angle (rad)')
title('tool1 input grasp angle')

fignum = fignum + 1;
figure(fignum)
scatter(mimic.t,mimic.tool{2}.master.grasp)
xlabel('time (s)')
ylabel('grap angle (rad)')
title('tool2 input grasp angle')

%% master vs slave motion
fontS = 14;

fignum = fignum + 1;
figure(fignum)
for i = 1:mimic.steps
   if(mimic.is_touch(i) == 1)
      if( mimic.tool{1}.slave.jaw_interact(i) == 17 || mimic.tool{1}.slave.jaw_interact(i) == 5) %rings or jacks
          h1 = quiver3(mimic.tool{1}.slave.X(i),mimic.tool{1}.slave.Y(i),mimic.tool{1}.slave.Z(i),mimic.tool{1}.slave.EA(i,1),mimic.tool{1}.slave.EA(i,2),mimic.tool{1}.slave.EA(i,3),0.05,'k');
          h2 = scatter3(mimic.tool{1}.master.X(i),mimic.tool{1}.master.Y(i),mimic.tool{1}.master.Z(i),'bo');
      end
   end
   hold on
end
hold off
grid on
title('Master vs slave tool motion')
xlabel('x (cm)')
ylabel('y (cm)')
zlabel('z (cm)')
legend([h1, h2],'Tool1 grasping (right)','Tool1 master')


%tool 2
fignum = fignum + 1;
figure(fignum)
for i = 1:mimic.steps
   if(mimic.is_touch(i) == 1)
      if( mimic.tool{2}.slave.jaw_interact(i) == 17 || mimic.tool{2}.slave.jaw_interact(i) == 5)
          h1 = quiver3(mimic.tool{2}.slave.X(i),mimic.tool{2}.slave.Y(i),mimic.tool{2}.slave.Z(i),mimic.tool{2}.slave.EA(i,1),mimic.tool{2}.slave.EA(i,2),mimic.tool{2}.slave.EA(i,3),0.05,'r');
          h2 = scatter3(mimic.tool{2}.master.X(i),mimic.tool{2}.master.Y(i),mimic.tool{2}.master.Z(i),'bo');
      end 
   end
   hold on
end
hold off
grid on
title('master vs slave tool motion')
xlabel('x (cm)')
ylabel('y (cm)')
zlabel('z (cm)')
legend([h1, h2],'Tool2 grasping (left)','Tool2 master')



%% Total motion
fontS = 14;

fignum = fignum + 1;
figure(fignum)
% segment motion

for i = 1:mimic.steps
   if(mimic.is_touch(i) == 1)
      if( mimic.tool{1}.slave.jaw_interact(i) == 17 || mimic.tool{1}.slave.jaw_interact(i) == 5)
          h1 = quiver3(mimic.tool{1}.slave.X(i),mimic.tool{1}.slave.Y(i),mimic.tool{1}.slave.Z(i),mimic.tool{1}.slave.EA(i,1),mimic.tool{1}.slave.EA(i,2),mimic.tool{1}.slave.EA(i,3),0.1,'r');
      end 
      hold on
      if( mimic.tool{2}.slave.jaw_interact(i) == 17 || mimic.tool{2}.slave.jaw_interact(i) == 5)
          h2 = quiver3(mimic.tool{2}.slave.X(i),mimic.tool{2}.slave.Y(i),mimic.tool{2}.slave.Z(i),mimic.tool{2}.slave.EA(i,1),mimic.tool{2}.slave.EA(i,2),mimic.tool{2}.slave.EA(i,3),0.1,'k');
      end
   end
   hold on
end
hold off
grid on
title('Segmented Tool Trajectory','FontSize',fontS)
xlabel('X (cm)','FontSize',fontS)
ylabel('Y (cm)','FontSize',fontS)
zlabel('Z (cm)','FontSize',fontS)
h_legend=legend([h1, h2],'Tool1 grasping (right)','Tool2 grasping (left)')
set(h_legend,'FontSize',12);

%% Camera clutching

fignum = fignum + 1;
figure(fignum)
% segment motion

for i = 1:mimic.steps
   if(mimic.pedals(i) == 6)
        h1 = scatter3(mimic.tool{1}.master.X(i),mimic.tool{1}.master.Y(i),mimic.tool{1}.master.Z(i),'ro');   
        hold on
        h2 = scatter3(mimic.tool{2}.master.X(i),mimic.tool{2}.master.Y(i),mimic.tool{2}.master.Z(i),'ko');
   end
   hold on
end

hold off
grid on
title('camera clutch - master motion')
xlabel('x (cm)')
ylabel('y (cm)')
zlabel('z (cm)')
legend([h1, h2],'Input1 clutch (right)','input2 clutch (left)')


%% plot velocities

fignum = fignum + 1;
figure(fignum)
plot(mimic.t,mimic.tool{1}.slave.dX,'rx')
title('tool 1 dx')

fignum = fignum + 1;
figure(fignum)
plot(mimic.t,mimic.tool{1}.slave.dY,'rx')
title('tool 1 dy')

fignum = fignum + 1;
figure(fignum)
plot(mimic.t,mimic.tool{1}.slave.dZ,'rx')
title('tool 1 dz')

%%

fignum = fignum + 1;
figure(fignum)
plot(mimic.t,mimic.tool{1}.slave.EA(:,1),'rx')
title('tool 1 roll')

fignum = fignum + 1;
figure(fignum)
plot(mimic.t,mimic.tool{1}.slave.EA(:,2),'rx')
title('tool 1 pitch')

fignum = fignum + 1;
figure(fignum)
plot(mimic.t,mimic.tool{1}.slave.EA(:,3),'rx')
title('tool 1 yaw')

%%
fignum = fignum + 1;
figure(fignum)
plot(1:length(mimic.deltaT),mimic.deltaT,'r-')
title('dt time fluxuation')