
close all

%% importing raw data
clear all

%plot settings
fontS = 12;

endData = 21000;
countsPerRev = 500; %encoder (astrini thesis)

filename = {'bladder4Hz1.lvm';'bladder4Hz2.lvm';'gallbladder4Hz1.lvm';'gallbladder4Hz2.lvm';'liver4Hz1.lvm';'liver4Hz2.lvm'};
enumTypes = {'Bladder'; 'Gallbladder'; 'Liver'};% 'smallbowel' };
enumAbbrev = {'B'; 'GB'; 'L'};% 'SB' };
tissuetypeIn = [1,1,2,2,3,3];

numFiles = length(filename);

headerlinesIn = 24;

for kk = 1:numFiles
    rawd = [];
    
    % data = [ t stress motorCommandActual jaw1contact jaw2contact encoder ] ;
    [tt ,rawd(:,1),rawd(:,2),rawd(:,3),rawd(:,4),rawd(:,5)]  =  textread(char(filename(kk)),'%f %f %f %f %f %f' ,'headerlines', headerlinesIn );

    %angle of grasper
    angle = rawd(:,5);
    angle = (angle / countsPerRev) * 2*pi; %convert to radians.
    
    %time
    time = tt;
    deltaT = mean(diff(time));

    FullData{kk}.t = time;
    FullData{kk}.tissuetype = tissuetypeIn(kk);
    FullData{kk}.ind = 1:1:length(time);
    FullData{kk}.stress = rawd(:,1); %newtons?
    FullData{kk}.stressdot = Calculate_velocity( FullData{kk}.stress, deltaT, 'holobrodko');
    FullData{kk}.stressdotdot = Calculate_velocity( FullData{kk}.stressdot, deltaT, 'holobrodko');
    FullData{kk}.angle = angle; %radians
    FullData{kk}.angledot = Calculate_velocity( FullData{kk}.angle, deltaT, 'holobrodko');  % FAR better than diff USE IT
    FullData{kk}.angledotdot = Calculate_velocity( FullData{kk}.angledot, deltaT, 'holobrodko');

    %In order to see directions
    plotlength = length(FullData{kk}.angle); %8000;
    figure(kk)
    
    subplot(3,1,1)
    plot(FullData{kk}.angle(1:plotlength))
    title('angle')
    
    subplot(3,1,2)
    plot(FullData{kk}.stress(1:plotlength))
    title('stress')

    subplot(3,1,3)
    plot(FullData{kk}.angledot(1:plotlength))
    title('angledot')
    
    str = sprintf('Raw Data: Grasp #: %i, tissue: %s',kk, enumTypes{ FullData{kk}.tissuetype } );
    suptitle(str);
    
end

%% GO BACKWARD IN TIME!

defgraspstep = 0.2;

% SIGNS ARE CORRECT FOR RAW
stressThreshold = 1.2; %Newtons, Critical threshold at which to start counting a grasp. As small as possible!
stressThreshHigh = 3; %Newtons, Critical threshold at which to start counting a grasp. As small as possible!
VelThreshold = 1; %mm/s, minimum positive closing speed at which to consider grasp as occuring
MinGraspInterval = 30; %Minimum number of time steps between grasps
MinGraspStep = defgraspstep; %Minimum namount of time a grasp should last
MinAngleDisp = 1;

%Go through every file
for kk = 1:numFiles

    %start grasp at end
    GraspStart_i = length(FullData{kk}.t);
    GraspEnd_i = length(FullData{kk}.t);
    ss = GraspEnd_i;
    ee = GraspEnd_i;
    
    %count up total grasps
    nGrasp = 0;
    
    %while stuff is available
    stillEnd = true;
    stillStart = true;
    while(stillEnd)
        
        %back up from last grasp start to find the end of the previous
        ee = ee - 1;
        if(ee < 1)
            stillEnd = false;
            break;
        end
        
        %reset while
        stillStart = true;
        
        %change min grasp step for last grasp
        if(nGrasp == 0)
           MinGraspStep = 0;
        else
           MinGraspStep = defgraspstep;
        end
        
        % Find the end of the grasp, top of the peak
        if (FullData{kk}.t(ee) < FullData{kk}.t(max(GraspStart_i - MinGraspInterval,1)) ) && (FullData{kk}.stress(ee) > stressThreshHigh) && ( FullData{kk}.angledot(ee) > VelThreshold )
            GraspEnd_i = ee;
            ss = GraspEnd_i;
            
            %Now find the start of that grasp
            while(stillStart) 
                
                %back up from grasp end to find start
                ss = ss - 1;
                if(ss < 1)
                    stillStart = false;
                    stillEnd = false;
                    break;
                end
                
                %Find index of the start of that grasp
                if (FullData{kk}.t(ss) < FullData{kk}.t(GraspEnd_i) - MinGraspStep) && (FullData{kk}.stress(ss) < stressThreshold) && (FullData{kk}.angle(GraspEnd_i) - FullData{kk}.angle(ss) > MinAngleDisp)
                    
                    GraspStart_i = ss;
                    
                    nGrasp = nGrasp + 1;
                    segDataTemp{kk}.GraspData{nGrasp}.t = FullData{kk}.t(GraspStart_i:GraspEnd_i);
                    segDataTemp{kk}.GraspData{nGrasp}.stress = FullData{kk}.stress(GraspStart_i:GraspEnd_i);
                    segDataTemp{kk}.GraspData{nGrasp}.stressdot = FullData{kk}.stressdot(GraspStart_i:GraspEnd_i);
                    segDataTemp{kk}.GraspData{nGrasp}.stressdotdot = FullData{kk}.stressdotdot(GraspStart_i:GraspEnd_i);
                    segDataTemp{kk}.GraspData{nGrasp}.angle = FullData{kk}.angle(GraspStart_i:GraspEnd_i);
                    segDataTemp{kk}.GraspData{nGrasp}.angledot = FullData{kk}.angledot(GraspStart_i:GraspEnd_i);
                    segDataTemp{kk}.GraspData{nGrasp}.angledotdot = FullData{kk}.angledotdot(GraspStart_i:GraspEnd_i);
                    segDataTemp{kk}.GraspData{nGrasp}.sampleNum = length( segDataTemp{kk}.GraspData{nGrasp}.t );
                    segDataTemp{kk}.GraspData{nGrasp}.tissuetype = FullData{kk}.tissuetype;
                    
                    fprintf('end: %f, start: %f \n',ee,ss);
                    
                    ee = GraspStart_i;
                    stillStart = false;
                end
                
            end

        end
        
    end
end


%Now reorder so grasps go in correct direction vs index
for kk = 1:numFiles
    %%%%%%%%%%%%%%%%%%%%%%%%%% segment grasps %%%%%%%%%%%%%%%%%%%%%%%
    segData{kk}.t = [];
    segData{kk}.stress = [];
    segData{kk}.stressdot = [];
    segData{kk}.stressdotdot = [];
    segData{kk}.angle = [];
    segData{kk}.angledot = [];
    segData{kk}.angledotdot = [];

    
    indnew = 1;
    %flip the grasp order so linearly increases in time
    for jj = length(segDataTemp{kk}.GraspData):-1:1
        segData{kk}.GraspData{indnew} = segDataTemp{kk}.GraspData{jj};
        
        %all grasps for plots
        segData{kk}.t = [segData{kk}.t; segData{kk}.GraspData{indnew}.t];
        segData{kk}.stress = [segData{kk}.stress; segData{kk}.GraspData{indnew}.stress];
        segData{kk}.stressdot = [segData{kk}.stressdot; segData{kk}.GraspData{indnew}.stressdot];
        segData{kk}.stressdotdot = [segData{kk}.stressdotdot; segData{kk}.GraspData{indnew}.stressdotdot];
        segData{kk}.angle = [segData{kk}.angle; segData{kk}.GraspData{indnew}.angle];
        segData{kk}.angledot = [segData{kk}.angledot; segData{kk}.GraspData{indnew}.angledot];
        segData{kk}.angledotdot = [segData{kk}.angledotdot; segData{kk}.GraspData{indnew}.angledotdot];
        segData{kk}.tissuetype = FullData{kk}.tissuetype;
        
        indnew = indnew + 1;
    end
    
    figure(kk)
    subplot(3,1,1)
    plot(FullData{kk}.t,FullData{kk}.angle)
    hold on
    scatter(segData{kk}.t,segData{kk}.angle,'.r')
    hold off
    xlabel('Time (s)')
    ylabel('angle (rad)')
    
    subplot(3,1,2)
    plot(FullData{kk}.t,FullData{kk}.angledot)
    hold on
    scatter(segData{kk}.t,segData{kk}.angledot,'.r')
    hold off
    xlabel('Time (s)')
    ylabel('angledot (rad)')
    
    subplot(3,1,3)
    plot(FullData{kk}.t,FullData{kk}.stress)
    hold on
    scatter(segData{kk}.t,segData{kk}.stress,'.r')
    hold off
    xlabel('Time (s)')
    ylabel('stress (N)')
    
    str = sprintf('Segmented Data: Grasp #: %i, tissue: %s',kk, enumTypes{ FullData{kk}.tissuetype } );
    suptitle(str);
    
end

return

%% subplot angle, angledot, angledotdot

for kk = 1:numFiles
    
    figure(kk)

    subplot(3,1,1)
    plot(FullData{kk}.t,FullData{kk}.angle)
    hold on
    scatter(segData{kk}.t, segData{kk}.angle,'.r')
    hold off
    
    %strtitle = strcat('Grasp Angle. Tissue = ', enumTypes(FullData{kk}.tissuetype ));
    title('Grasp Angle')
    xlabel('Time (s)')
    ylabel('Angle (rad)')

    subplot(3,1,2)
    plot(FullData{kk}.t,FullData{kk}.angledot)
    hold on
    scatter(segData{kk}.t, segData{kk}.angledot,'.r')
    hold off
    
    %strtitle = strcat('Angledot. Tissue = ', enumTypes(FullData{kk}.tissuetype ));
    title('Angledot')
    xlabel('Time (s)')
    ylabel('Angledot (rad/s)')

    subplot(3,1,3)
    plot(FullData{kk}.t,FullData{kk}.angledotdot)
    hold on
    scatter(segData{kk}.t, segData{kk}.angledotdot,'.r')
    hold off
    
    %strtitle = strcat('Angledotdot. Tissue = ', enumTypes(FullData{kk}.tissuetype ));
    title('Angledotdot')
    xlabel('Time (s)')
    ylabel('Angledotdot (rad/s2)')
    
    str = sprintf('Segmented Grasp Angle: Grasp #: %i, tissue: %s',kk, enumTypes{ FullData{kk}.tissuetype } );
    suptitle(str);

end

%% Phase Portraits
[valind, typeIndex] = unique(tissuetypeIn); %first unique number index

figure
CM = { [1,0,0] , [0,1,0], [0,0,1], [1,1,0], [0,1,1], [1,0,1], [0,0,0], [.5,1,1] };

for kk = 1:numFiles
    h(tissuetypeIn(kk)) = quiver(segData{kk}.angle,segData{kk}.angledot,segData{kk}.angledot,segData{kk}.angledotdot,'color',CM{FullData{kk}.tissuetype});
    
    %legend_name{kk} = enumTypes{ FullData{kk}.tissuetype };
    hold on
end
hold off

title('Phase Portrait (angle, angledot)')
xlabel('Angle (rad)')
ylabel('Angledot (rad/s)')
legend([h(1),h(2),h(3),h(4)],enumTypes{:})

figure

for kk = 1:numFiles
    g(tissuetypeIn(kk)) = quiver(segData{kk}.stress,segData{kk}.stressdot,segData{kk}.stressdot,segData{kk}.stressdotdot,'color',CM{FullData{kk}.tissuetype});

    %legend_name{kk} = enumTypes{ FullData{kk}.tissuetype };
    hold on
end
hold off

title('Phase Portrait (stress, stressdot)')
xlabel('stress (N)')
ylabel('stressdot (N/s)')
legend([g(1),g(2),g(3),g(4)],enumTypes{:})


figure

for kk = 1:numFiles
    f(tissuetypeIn(kk)) = quiver(segData{kk}.stress,segData{kk}.angledot,segData{kk}.stressdot,segData{kk}.angledotdot,'color',CM{FullData{kk}.tissuetype});
    
    %legend_name{kk} = enumTypes{ FullData{kk}.tissuetype };
    hold on
end
hold off

title('Phase Portrait (stress, angledot)')
xlabel('stress (N)')
ylabel('angledot (rad/s)')
legend([f(1),f(2),f(3),f(4)],enumTypes{:})

figure

for kk = 1:numFiles
    f(kk) = quiver(segData{kk}.stress,segData{kk}.angledot,segData{kk}.stressdot,segData{kk}.angledotdot,'color',CM{kk});
    
    %legend_name{kk} = enumTypes{ FullData{kk}.tissuetype };
    hold on
end
hold off

title('Phase Portrait files (stress, angledot)')
xlabel('stress (N)')
ylabel('angledot (rad/s)')
legend([f(1),f(2),f(3),f(4),f(5),f(6),f(7),f(8)],filename{:})


%% 3D Phase Portraits
[valind, typeIndex] = unique(tissuetypeIn); %first unique number index

figure
CM = { [1,0,0] , [0,1,0], [0,0,1], [1,1,0], [0,1,1], [1,0,1], [0,0,0], [.5,1,1] };


sss = 15;

for kk = 1:numFiles
    h(tissuetypeIn(kk)) = scatter3(segData{kk}.angle,segData{kk}.angledot,segData{kk}.angledotdot,sss,CM{FullData{kk}.tissuetype});
    
    %legend_name{kk} = enumTypes{ FullData{kk}.tissuetype };
    hold on
end
hold off

title('Phase Portrait')
xlabel('Angle (rad)')
ylabel('Angledot (rad/s)')
zlabel('Angledotdot (rad/s^2)')
legend([h(1),h(2),h(3)],enumTypes{:})


%% combine segments

%set to empty
for jj = 1:max(tissuetypeIn)
    combineSeg{jj}.t = [];
    combineSeg{jj}.stress = [];
    combineSeg{jj}.stressdot = [];
    combineSeg{jj}.stressdotdot = [];
    combineSeg{jj}.angle = [];
    combineSeg{jj}.angledot = [];
    combineSeg{jj}.angledotdot = [];
    ombineSeg{jj}.D = [];
    combineSeg{jj}.U = [];
    combineSeg{jj}.Phi = [];
end

% combined segmented grasps
for kk = 1:numFiles
    combineSeg{FullData{kk}.tissuetype}.t = [combineSeg{FullData{kk}.tissuetype}.t ; segData{kk}.t ];
    combineSeg{FullData{kk}.tissuetype}.stress = [combineSeg{FullData{kk}.tissuetype}.stress ; segData{kk}.stress ];
    combineSeg{FullData{kk}.tissuetype}.stressdot = [combineSeg{FullData{kk}.tissuetype}.stressdot ; segData{kk}.stressdot ];
    combineSeg{FullData{kk}.tissuetype}.stressdotdot = [combineSeg{FullData{kk}.tissuetype}.stressdotdot ; segData{kk}.stressdotdot ];
    combineSeg{FullData{kk}.tissuetype}.angle = [combineSeg{FullData{kk}.tissuetype}.angle ; segData{kk}.angle ];
    combineSeg{FullData{kk}.tissuetype}.angledot = [combineSeg{FullData{kk}.tissuetype}.angledot ; segData{kk}.angledot ];
    combineSeg{FullData{kk}.tissuetype}.angledotdot = [combineSeg{FullData{kk}.tissuetype}.angledotdot ; segData{kk}.angledotdot ];
    
end


lambda = [0:0.1:1];

%Create D and U matrix for RLS
for jj = 1:max(tissuetypeIn)
    combineSeg{jj}.state = [combineSeg{jj}.angledotdot, combineSeg{jj}.angledot, combineSeg{jj}.angle, combineSeg{jj}.angle.^2, combineSeg{jj}.angle.^3, ones(length(combineSeg{jj}.angle) , 1) ];
    combineSeg{jj}.input = [combineSeg{jj}.stress ];
    combineSeg{jj}.Phi = inv(combineSeg{jj}.state'*combineSeg{jj}.state) * combineSeg{jj}.state' * combineSeg{jj}.input;
end


%Discriminant least squares
for kk = 1:length(lambda)
    DLS_Parameters{kk} = DLS_TrainGeneral(combineSeg, ones(1,max(tissuetypeIn))*lambda(kk), 0);
end

%should evaluate to:
phi1 = [-0.000621406370117322;0.0976101250736223;-0.388256674969638;0.131870191736068;0.0373350034929734];
phi2 = [-0.000826176494778239;0.101152437533438;0.503038080455392;-0.429175896084507;0.129236022201036];
phi3 = [-0.00174977066715285;0.131433006396432;-0.246563945926546;0.0476950832392320;0.0553277583906477];
phi4 = [-0.000617185498974871;0.0233326316387804;0.679253921276928;-0.238591813349058;0.0703108547058125];


%% Online


lambda = [0:0.01:0.1];
yes_all_lambdas = 2;
classifyTotal = [];

%zero out classificiation struct
if(yes_all_lambdas == 1)
    for kk = 1:numFiles
        for tt = 1:length(lambda)
            classify{segData{kk}.tissuetype}.lamba{tt} = [];
        end
    end
end
if(yes_all_lambdas == 2)
    for kk = 1:numFiles
        for tt1 = 1:length(lambda) 
            for tt2 = 1:length(lambda) 
                for tt3 = 1:length(lambda) 
                    %for tt4 = 1:length(lambda) 
                        classify{segData{kk}.tissuetype}.l1{tt1}.l2{tt2}.l3{tt3} = [];
                    %end
                end
            end
        end
    end
end

ind = 1;

tic

%loop through all files to get the leave one out
for kk = 1:numFiles
    for jj = 1:length( segData{kk}.GraspData )
        %Clear mats
        for pp = 1:max(tissuetypeIn)
            trainSeg{pp}.state = [];
            trainSeg{pp}.input = [];
        end
        online.state = [];
        online.input = [];
        %loop through actually all files
        for ii = 1:numFiles
            for hh = 1:length( segData{ii}.GraspData )
                %Populate training and online data sets
                
                if(ii == kk)
                    % if current file is slated for leave-one-out
                    if(hh == jj)
                        %if current grasp is slated for leave-one-out 
                        %(add it to online)
                        online.state = [segData{ii}.GraspData{hh}.angledotdot, segData{ii}.GraspData{hh}.angledot, segData{ii}.GraspData{hh}.angle, segData{ii}.GraspData{hh}.angle.^2, segData{ii}.GraspData{hh}.angle.^3, ones(length(segData{ii}.GraspData{hh}.angle) , 1) ];
                        online.input = [segData{ii}.GraspData{hh}.stress ];
                    else
                        %current grasp is just for training
                        trainSeg{FullData{ii}.tissuetype}.state = [trainSeg{FullData{ii}.tissuetype}.state; segData{ii}.GraspData{hh}.angledotdot, segData{ii}.GraspData{hh}.angledot, segData{ii}.GraspData{hh}.angle, segData{ii}.GraspData{hh}.angle.^2, segData{ii}.GraspData{hh}.angle.^3, ones(length(segData{ii}.GraspData{hh}.angle) , 1) ];
                        trainSeg{FullData{ii}.tissuetype}.input = [trainSeg{FullData{ii}.tissuetype}.input; segData{ii}.GraspData{hh}.stress ];
                    end
                else
                    % if current file is just for training
                    trainSeg{FullData{ii}.tissuetype}.state = [trainSeg{FullData{ii}.tissuetype}.state; segData{ii}.GraspData{hh}.angledotdot, segData{ii}.GraspData{hh}.angledot, segData{ii}.GraspData{hh}.angle, segData{ii}.GraspData{hh}.angle.^2, segData{ii}.GraspData{hh}.angle.^3, ones(length(segData{ii}.GraspData{hh}.angle) , 1) ];
                    trainSeg{FullData{ii}.tissuetype}.input = [trainSeg{FullData{ii}.tissuetype}.input; segData{ii}.GraspData{hh}.stress ];
                end

                
            end
        end
        
        %Discriminant least squares (All the same Lambdas)
        if(yes_all_lambdas == 1)
            for tt = 1:length(lambda) %for now all lambdas are the same
                tt
                DLS_LOO{tt} = DLS_TrainGeneral(trainSeg, ones(1,max(tissuetypeIn))*lambda(tt), 0);

                %Now Test Online classification
                for gg = 1:max(tissuetypeIn)
                    error_temp(gg) = sum(abs(online.input - online.state*DLS_LOO{tt}(:,gg) ));
                end

                %Get the minimum sum error
                [minner,idx] = min(error_temp);
                %save actual class and estimated class
                classify{segData{kk}.tissuetype}.lamba{tt} = [classify{segData{kk}.tissuetype}.lamba{tt}; segData{kk}.tissuetype, idx];

                classifyTotal(ind,:) = [segData{kk}.tissuetype, idx,lambda(tt)];
                ind = ind + 1;
            end
        end
        
        %Discriminant least squares (pairwise lambda combos
        if(yes_all_lambdas == 2)
            for tt1 = 1:length(lambda) 
                for tt2 = 1:length(lambda) 
                    for tt3 = 1:length(lambda) 
                        %for tt4 = 1:length(lambda) 

                            %DLS_LOO.l1{tt1}.l2{tt2}.l3{tt3}.l4{tt4} = DLS_TrainGeneral(trainSeg, [lambda(tt1),lambda(tt2),lambda(tt3),lambda(tt4)], 0);
                            DLS_LOO.l1{tt1}.l2{tt2}.l3{tt3} = DLS_TrainGeneral(trainSeg, [lambda(tt1),lambda(tt2),lambda(tt3)], 0);

                            %Now Test Online classification
                            for gg = 1:max(tissuetypeIn)
                                error_temp(gg) = sum(abs(online.input - online.state*DLS_LOO.l1{tt1}.l2{tt2}.l3{tt3}(:,gg) ));
                            end

                            %Get the minimum sum error
                            [minner,idx] = min(error_temp);
                            %save actual class and estimated class
                            classify{segData{kk}.tissuetype}.l1{tt1}.l2{tt2}.l3{tt3} = [classify{segData{kk}.tissuetype}.l1{tt1}.l2{tt2}.l3{tt3}; segData{kk}.tissuetype, idx];

                            classifyTotal(ind,:) = [segData{kk}.tissuetype, idx, lambda(tt1),lambda(tt2),lambda(tt3)];
                            ind = ind + 1;
                        %end
                    end
                end
            end
        end
        
    end
end

toc

%Get classification precentage
correct = (classifyTotal(:,1) == classifyTotal(:,2) );
disp('classification percentage: ')
corr_percent = sum(correct)/length(correct)

%% Find best lambdas

best_percent = 0;
best_lambda = [];
corr_all = [];


for tt1 = 1:length(lambda) 
    for tt2 = 1:length(lambda) 
        for tt3 = 1:length(lambda) 
            %for tt4 = 1:length(lambda)
                for ii = 1:max(tissuetypeIn)
                    correct_all = classify{ii}.l1{tt1}.l2{tt2}.l3{tt3}(:,1) == classify{ii}.l1{tt1}.l2{tt2}.l3{tt3}(:,2);
                    corr_percent(ii) = sum(correct_all)/length(correct_all);
                end
                
                
                
                best_new = mean(corr_percent);
                if(best_new > best_percent)
                    best_percent = best_new;
                    best_lambda = [lambda(tt1), lambda(tt2), lambda(tt3)];
                end
                
                corr_all = [corr_all; best_new, corr_percent, lambda(tt1), lambda(tt2), lambda(tt3),tt1,tt2,tt3];
            %end
        end
    end
end

sort_col = 1;
[Y_all,I_all]=sort(corr_all(:,sort_col));
corr_sorted = corr_all(I_all,:); 

disp('best lambdas:')

lambdars = corr_sorted(end,5:7)

disp('best phis:')

Phi_Best = DLS_TrainGeneral(combineSeg, [corr_sorted(end,5),corr_sorted(end,6),corr_sorted(end,7)], 0)


%% heat map

corr_heat = corr_sorted(:,[1,5,6,7]);

G_types = [1,2,3];

FIG_W = 3.5;     % Width of actual figure  
FIG_H = 2;     % Height of actual figure
FIG_UNITS = 'inches'; % units for W&H
FIG_RES = 400; % figure resolution in dpi

for ii = 1:max(tissuetypeIn)
    for jj = 1:max(tissuetypeIn)
        in_types = [ii,jj];
        if(ii ~= jj)
            plotter = corr_heat(corr_heat(:,jj+1) == 0,:)
            cross = setxor(G_types, in_types); %find whats left
            
            %get xyz for plot
            x = plotter(:,ii+1);
            y = plotter(:,cross+1);
            z = plotter(:,1)*100;
            
            %trick to get heatmap from xyz
            xlin = linspace(min(x),max(x),33);
            ylin = linspace(min(y),max(y),33);
            [X,Y] = meshgrid(xlin,ylin);
            f = scatteredInterpolant(x,y,z);
            Z = f(X,Y);
            
            figure
            %mesh(X,Y,Z) %interpolated
            contourf(X,Y,Z)
            axis tight; hold on
            colormap(cool);
            c = colorbar;
            c.Label.String = 'Accuracy (%)';
            str = sprintf('Classification vs \\lambda, %s and %s', enumTypes{ ii } , enumTypes{ cross });
            title(str,'FontSize',10,'FontName','Times')
            strx = sprintf('\\lambda_{%s}', enumAbbrev{ ii });
            stry = sprintf('\\lambda_{%s}', enumAbbrev{ cross });
            xlabel(strx,'FontSize',10,'FontName','Times')
            ylabel(stry,'FontSize',10,'FontName','Times')
            
            
            % set it's W and H w/o messing up the position on the screen
            set(gcf,'PaperPositionMode','auto', 'units', FIG_UNITS)
            FIG_SZ = get(gcf, 'position');
            FIG_SZ(3:end) = [FIG_W FIG_H];
            set(gcf, 'position', FIG_SZ);

        end
        
    end
end

% figure
% scatter3(corr_sorted(:,5),corr_sorted(:,6),corr_sorted(:,7),8,corr_sorted(:,1));


%% convergence rate

ind = 1;

%Create D and U matrix for RLS
for jj = 1:max(tissuetypeIn)
    combineSegC{jj}.state = [combineSeg{jj}.angledotdot, combineSeg{jj}.angledot, combineSeg{jj}.angle, combineSeg{jj}.angle.^2, combineSeg{jj}.angle.^3, ones(length(combineSeg{jj}.angle) , 1) ];
    combineSegC{jj}.input = [combineSeg{jj}.stress ];
end
        
DLS_ALPHA = DLS_TrainGeneral(combineSegC, lambdars, 0);

for cc = 1:max(tissuetypeIn)
    ngrasp(cc) = 1;
end

%Now Test Online classification
for ii = 1:numFiles
    for hh = 1:length( segData{ii}.GraspData )
        online.state = [segData{ii}.GraspData{hh}.angledotdot, segData{ii}.GraspData{hh}.angledot, segData{ii}.GraspData{hh}.angle, segData{ii}.GraspData{hh}.angle.^2, segData{ii}.GraspData{hh}.angle.^3, ones(length(segData{ii}.GraspData{hh}.angle) , 1) ];
        online.input = [segData{ii}.GraspData{hh}.stress ];
        errorsums{segData{ii}.tissuetype}.grasp{ngrasp(segData{ii}.tissuetype)} = [];
        for gg = 1:max(tissuetypeIn)
            
            errorsums{segData{ii}.tissuetype}.grasp{ngrasp(segData{ii}.tissuetype)} = [errorsums{segData{ii}.tissuetype}.grasp{ngrasp(segData{ii}.tissuetype)}, cumsum(abs(online.input - online.state*DLS_ALPHA(:,gg) )) ];
        end
        ngrasp(segData{ii}.tissuetype) = ngrasp(segData{ii}.tissuetype) + 1;
    end
end

for cc = 1:max(tissuetypeIn)
    corr_time = [];
    for ii = 1:ngrasp(cc)-1
        [vallrs,idxmin] = min(errorsums{cc}.grasp{ii},[],2)
        correct_all = idxmin == cc;
        corr_time = [corr_time , sum(correct_all)/length(correct_all)];
    end
    corr_time_mean(cc) = mean(corr_time);
end

corr_time_mean

alpha_thresh = 0.5;

ALPHA_SCALE(1) = 10;
ALPHA_SCALE(2) = 10;
ALPHA_SCALE(3) = 10;


FIG_W = 3.5;     % Width of actual figure  
FIG_H = 2;     % Height of actual figure
FIG_UNITS = 'inches'; % units for W&H
FIG_RES = 400; % figure resolution in dpi

%Plot alphas
figure
cc = 1;
for ii = 1:ngrasp(cc)-1
    diffe = abs( ( (errorsums{cc}.grasp{ii}(:,2) .* errorsums{cc}.grasp{ii}(:,3)) - errorsums{cc}.grasp{ii}(:,cc) )  )./ ( errorsums{cc}.grasp{ii}(:,2) .* errorsums{cc}.grasp{ii}(:,3) + ALPHA_SCALE(cc) );
    dt = 0.001;
    h1 = plot([1:length(diffe)]*dt,diffe,'r-');
    idxc{cc}(ii) = find(diffe > alpha_thresh, 1, 'first');
    mean_error{cc}(ii) = mean(errorsums{cc}.grasp{ii}(:,cc));
    mean_diff{cc}(ii) = mean(abs(errorsums{cc}.grasp{ii}(:,cc) - errorsums{cc}.grasp{ii}(:,2) - errorsums{cc}.grasp{ii}(:,3)));
    hold on
end

cc = 2;
for ii = 1:ngrasp(cc)-1
    diffe = abs( ( (errorsums{cc}.grasp{ii}(:,1) .* errorsums{cc}.grasp{ii}(:,3)) - errorsums{cc}.grasp{ii}(:,cc) )  )./ ( errorsums{cc}.grasp{ii}(:,1) .* errorsums{cc}.grasp{ii}(:,3) + ALPHA_SCALE(cc) );
    dt = 0.001;
    h2 = plot([1:length(diffe)]*dt,diffe,'b-');
    idxc{cc}(ii)  = find(diffe > alpha_thresh, 1, 'first');
    mean_error{cc}(ii)  = mean(errorsums{cc}.grasp{ii}(:,cc));
    mean_diff{cc}(ii)  = mean(abs(errorsums{cc}.grasp{ii}(:,cc) - errorsums{cc}.grasp{ii}(:,1) - errorsums{cc}.grasp{ii}(:,3)));
    hold on
end

cc = 3;
for ii = 1:ngrasp(cc)-1
    diffe = abs( ( (errorsums{cc}.grasp{ii}(:,1) .* errorsums{cc}.grasp{ii}(:,2)) - errorsums{cc}.grasp{ii}(:,cc) )  )./ ( errorsums{cc}.grasp{ii}(:,1) .* errorsums{cc}.grasp{ii}(:,2) + ALPHA_SCALE(cc) );
    dt = 0.001;
    h3 = plot([1:length(diffe)]*dt,diffe,'g-');
    idxc{cc}(ii)  = find(diffe > alpha_thresh, 1, 'first');
    mean_error{cc}(ii)  = mean(errorsums{cc}.grasp{ii}(:,cc));
    mean_diff{cc}(ii)  = mean(abs(errorsums{cc}.grasp{ii}(:,cc) - errorsums{cc}.grasp{ii}(:,1) - errorsums{cc}.grasp{ii}(:,2)));
    hold on
end

hold off
axis([0,0.11,0,1])
xlabel('Time (s)','FontSize',10,'FontName','Times')
ylabel('\alpha','FontSize',10,'FontName','Times')
legend([h1(1),h2(1),h3(1)],{enumTypes{1},enumTypes{2},enumTypes{3}},'FontSize',8,'Location','southeast','FontName','Times')

% set it's W and H w/o messing up the position on the screen
set(gcf,'PaperPositionMode','auto', 'units', FIG_UNITS)
FIG_SZ = get(gcf, 'position');
FIG_SZ(3:end) = [FIG_W FIG_H];
set(gcf, 'position', FIG_SZ);


% Save the figure to file
%print(gcf, 'C:\Users\MRDLAB\Documents\Rod Dockter\Surgical Robotics\mrdPublications\ICRA_Letters_Grasper\Images\alpha_combined', '-dpng', ['-r' num2str(FIG_RES)]) % Raster

mean(idxc{1})
mean(idxc{2})
mean(idxc{3})


std(idxc{1})
std(idxc{2})
std(idxc{3})

mean(mean_error{1})
mean(mean_error{2})
mean(mean_error{3})
% mean_error
% mean_diff

mean(mean_diff{1})
mean(mean_diff{2})
mean(mean_diff{3})


%% Plot classification versus lambda for each class


% colormap = {[1 0 0],[0 1 0], [0 0 1], [0 1 1], [1 0 1], [0 1/2 1] , [1/2 0 1]};
% 
% figure
% 
% %zero out classificiation struct
% for kk = 1:max(tissuetypeIn)
%     for tt = 1:length(lambda)
%         correct_temp = (classify{kk}.lamba{tt}(:,1) == classify{kk}.lamba{tt}(:,2) );
%         percent_temp(tt) = (sum(correct_temp)/length(correct_temp)) * 100;
%     end
%     
%     h(kk) = plot(lambda,percent_temp,'-x','Color', colormap{kk},'LineWidth',3);
%     hold on
%     strngr{kk} = sprintf('Tissue Type = %s',enumTypes{kk});
% end
% 
% hold off
% 
% axis([min(lambda) max(lambda) 0 100])
% title('Discriminant Dynamics Improved Classification','FontSize',fontS)
% xlabel('% \lambda_{C}','FontSize',fontS)
% ylabel('% Correct Classification','FontSize',fontS)
% h_legend=legend(h,strngr);
% set(h_legend,'FontSize',fontS);
