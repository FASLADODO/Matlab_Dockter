%filename = {'gallbladder4Hz1.lvm'; 'gallbladder4Hz2.lvm';'gallbladder025Hz1.lvm';'gallbladder025Hz2.lvm';'smallbowel4Hz1.lvm';'smallbowel4Hz2.lvm';'smallbowel025Hz1.lvm';'smallbowel025Hz2.lvm'};

clear all

filename = {'puttycalibHardTissue.txt'; 'puttycalibSoftTissue.txt'};
enumTypes = {'Hard'; 'Soft' };
%tissuetypeIn = [1,1,1,1,2,2,2,2];
tissuetypeIn = [1,2];

numFiles = length(filename);

headerlinesIn = 24;

for kk = 1:numFiles
    rawd = [];
    
    
    
    % data = <angle, angledot, angledotdot, load cell, pot read>
    [tt ,rawd(:,1),rawd(:,2),rawd(:,3),rawd(:,4)]  =  textread(char(filename(kk)),'%f,%d,%f,%f,%d' );

    deltaT = mean(diff(tt));
    
    FullData{kk}.t = tt;
    FullData{kk}.tissuetype = tissuetypeIn(kk);
    FullData{kk}.ind = 1:1:length(tt);
    
    FullData{kk}.angle = rawd(:,1); %radians
    FullData{kk}.force = rawd(:,4); %newtons?
    
%     FullData{kk}.angledot = Calculate_velocity( FullData{kk}.angle, deltaT, 'holobrodko');  % FAR better than diff USE IT
%     FullData{kk}.angledotdot = Calculate_velocity( FullData{kk}.angledot, deltaT, 'holobrodko');
%     
    FullData{kk}.angledot = rawd(:,2);
    FullData{kk}.angledotdot = rawd(:,3);
    
    %force always non-negative
    FullData{kk}.force(FullData{kk}.force<0)=0;
    
    %In order to see directions
    plotlength = length(FullData{kk}.t);
    figure(kk)
    
    subplot(4,1,1)
    plot(FullData{kk}.t,FullData{kk}.angle(1:plotlength))
    xlabel('time (ms)')
    ylabel('angle (encoder)')
    strr = sprintf('Tissue type = %s',enumTypes{kk});
    title(strr);

    subplot(4,1,2)
    plot(FullData{kk}.t,FullData{kk}.angledot(1:plotlength))
    xlabel('time (ms)')
    ylabel('angledot (encoder/ms)')
    
    subplot(4,1,3)
    plot(FullData{kk}.t,FullData{kk}.angledotdot(1:plotlength))
    xlabel('time (ms)')
    ylabel('angledotdot (encoder/ms)')
    
    
    subplot(4,1,4)
    plot(FullData{kk}.t,FullData{kk}.force(1:plotlength))
    xlabel('time (ms)')
    ylabel('force (load cell)')

end

%% Segment Data

% SIGNS ARE CORRECT FOR RAW
forceThreshold = 20000; %Newtons, Critical threshold at which to start counting a grasp. As small as possible!
VelThreshold = -0.4; %mm/s, minimum positive closing speed at which to consider grasp as occuring
VelBackThreshold = 0.001; %mm/s, minimum positive closing speed at which to consider grasp as occuring
MinGraspInterval = 5; %Minimum number of time steps passed required before saying a grasp has ended
MinGraspStep = 60; %Minimum number of time steps passed required before saying a new grasp has started
MinAngleDisp = 25;

totalGrasps = [];

for kk = 1:numFiles 

    %%%%%%%%%%%%%%%%%%%%%%%%%% segment grasps %%%%%%%%%%%%%%%%%%%%%%%
    segData{kk}.t = [];
    
    segData{kk}.angle = [];
    segData{kk}.angledot = [];
    segData{kk}.angledotdot = [];
    segData{kk}.force = [];

    GraspStart_i = 1;
    GraspEnd_i = 1;
    nGrasp = 0;

    while ~isempty(GraspStart_i)
        if(nGrasp == 0)
           MinGraspStep = 0;
        else
           MinGraspStep = 100;
        end
        
        %Find index of where grasp begins
        %fprintf('start first &, %f \n',kk);
        GraspStart_i = find( (FullData{kk}.t > FullData{kk}.t(GraspStart_i)) & (FullData{kk}.t > FullData{kk}.t(GraspEnd_i) + MinGraspStep)  & (FullData{kk}.force > forceThreshold) & (FullData{kk}.angledot < VelThreshold) , 1,'first');
        %fprintf('first &, %f \n',kk);
        %Found grasp, lets add data
        if ~isempty(GraspStart_i)

            %Find index of the end of that grasp
            GraspEnd_i = find( (FullData{kk}.t > FullData{kk}.t(min(GraspStart_i + MinGraspInterval, length(FullData{kk}.t) ) )) & (FullData{kk}.angledot > VelBackThreshold), 1, 'first');
            %fprintf('second &, %f \n',kk);
            %It is okay to have decreasing force while xdot is still positive
            %because it indicates deceleration of grasp.        
            if FullData{kk}.angle(GraspStart_i) - FullData{kk}.angle(GraspEnd_i)  > MinAngleDisp
                nGrasp = nGrasp + 1;
                segData{kk}.GraspData{nGrasp}.t = FullData{kk}.t(GraspStart_i:GraspEnd_i);
                segData{kk}.GraspData{nGrasp}.force = FullData{kk}.force(GraspStart_i:GraspEnd_i);
                segData{kk}.GraspData{nGrasp}.angle = FullData{kk}.angle(GraspStart_i:GraspEnd_i);
                segData{kk}.GraspData{nGrasp}.angledot = FullData{kk}.angledot(GraspStart_i:GraspEnd_i);
                segData{kk}.GraspData{nGrasp}.angledotdot = FullData{kk}.angledotdot(GraspStart_i:GraspEnd_i);
                segData{kk}.GraspData{nGrasp}.sampleNum = length( segData{kk}.GraspData{nGrasp}.t );
                %fprintf('first seg, %f \n',kk);
                %fprintf('T diff=, %f \n',GraspEnd_i-GraspStart_i);
                %Add to allGrasps for plotting
                segData{kk}.t = [segData{kk}.t; segData{kk}.GraspData{nGrasp}.t];
                segData{kk}.force = [segData{kk}.force; segData{kk}.GraspData{nGrasp}.force];
                segData{kk}.angle = [segData{kk}.angle; segData{kk}.GraspData{nGrasp}.angle];
                segData{kk}.angledot = [segData{kk}.angledot; segData{kk}.GraspData{nGrasp}.angledot];
                segData{kk}.angledotdot = [segData{kk}.angledotdot; segData{kk}.GraspData{nGrasp}.angledotdot];
                %fprintf('second seg, %f \n',kk);
            end

            %Update to begin looking for next grasp
            GraspStart_i = GraspEnd_i + 5;

        end
        
        totalGrasps(kk) = nGrasp;
    end

    figure(kk)
    subplot(3,1,1)
    plot(FullData{kk}.t,FullData{kk}.angle)
    hold on
    scatter(segData{kk}.t,segData{kk}.angle,'.r')
    %scatter(GraspData{6}.t, GraspData{6}.angle,'.r')
    hold off
    xlabel('Time (s)')
    ylabel('angle (rad)')
    
    strr = sprintf('Segmented Grasps. Type = %s',enumTypes{kk});
    title(strr);

    
    subplot(3,1,2)
    plot(FullData{kk}.t,FullData{kk}.angledot)
    hold on
    scatter(segData{kk}.t,segData{kk}.angledot,'.r')
    %scatter(GraspData{6}.t, GraspData{6}.angle,'.r')
    hold off
    xlabel('Time (s)')
    ylabel('angledot (rad)')
    
    subplot(3,1,3)
    plot(FullData{kk}.t,FullData{kk}.force)
    hold on
    scatter(segData{kk}.t,segData{kk}.force,'.r')
    %scatter(GraspData{6}.t, GraspData{6}.angle,'.r')
    hold off
    xlabel('Time (s)')
    ylabel('Force (raw)')
end


%% Plot state space

CM = { [1,0,0] , [0,1,0], [0,0,1], [1,1,0], [0,1,1], [1,0,1], [0,0,0], [1,1,1] };
scale_factor = 10;
fontS = 14;

figure(69)

for kk = 1:numFiles
    h(kk) = quiver(segData{kk}.angle,segData{kk}.angledot,segData{kk}.angledot*scale_factor,segData{kk}.angledotdot*scale_factor,'color',CM{kk},'AutoScale','off');
    adjust_quiver_arrowhead_size(h(kk), 0.04);
    hold on
end
hold off

title('Sample tissue grasp phase portrait','FontSize',fontS)
xlabel('angle','FontSize',fontS)
ylabel('angle dot','FontSize',fontS)
legend([h(1),h(2)],enumTypes{1},enumTypes{2})


%% Solve params

D1 = [segData{1}.angledotdot, segData{1}.angledot, segData{1}.angle, segData{1}.angle.^2, segData{1}.angle.^3];
U1 = segData{1}.force;
Phi1 = inv(D1' * D1) * D1' * U1;
Phi1_prime = lsqr(D1,U1);

D2 = [segData{2}.angledotdot, segData{2}.angledot, segData{2}.angle, segData{2}.angle.^2, segData{2}.angle.^3];
U2 = segData{2}.force;
Phi2 = inv(D2' * D2) * D2' * U2;
Phi2_prime = lsqr(D2,U2);


%%
classify = [];
est_class = 0;
cntr = 0;
errors{1} = [];

for kk = 1:numFiles
    for jj = 1:totalGrasps(kk)
        cntr = cntr + 1;
        
        Dx = [segData{kk}.GraspData{jj}.angledotdot,segData{kk}.GraspData{jj}.angledot,segData{kk}.GraspData{jj}.angle,segData{kk}.GraspData{jj}.angle.^2,segData{kk}.GraspData{jj}.angle.^3];
        
        e1 = abs(segData{kk}.GraspData{jj}.force - Dx*Phi1_prime);
        e2 = abs(segData{kk}.GraspData{jj}.force - Dx*Phi2_prime);
        if(sum(e1) < sum(e2))
            est_class = 1;
        else
            est_class = 2;
        end
        
        classify(cntr,:) = [kk,est_class,sum(e1),sum(e2)];
        errors{kk}.grasp{jj} = [e1,e2];
        error_sums{kk}.grasp{jj} = [cumsum(e1),cumsum(e2)];
        
        if(kk == 1)
            arr = cumsum(e1) < cumsum(e2);
            timgrasp = segData{kk}.GraspData{jj}.t - segData{kk}.GraspData{jj}.t(1);
            convergence{kk}.comp{jj} = [arr, timgrasp];
            idx = find(arr == 0, 1, 'last');
            if(~isempty(idx))
                convergence{kk}.time(jj) = timgrasp(idx+1);
            else
                convergence{kk}.time(jj) = 0;
            end
        else
            arr = cumsum(e2) < cumsum(e1);
            timgrasp = segData{kk}.GraspData{jj}.t - segData{kk}.GraspData{jj}.t(1);
            convergence{kk}.comp{jj} = [arr, timgrasp];
            idx = find(arr == 0, 1, 'last');
            if(~isempty(idx))
                convergence{kk}.time(jj) = timgrasp(idx+1);
            else
                convergence{kk}.time(jj) = 0;
            end
        end
    end
end

average_convergence_class_1 = mean(convergence{1}.time)
average_convergence_class_2 = mean(convergence{2}.time)

figure(11)
for jj = 1:totalGrasps(1)
    alphas = ( error_sums{1}.grasp{jj}(:,2) - error_sums{1}.grasp{jj}(:,1) ) ./ error_sums{1}.grasp{jj}(:,2);
    
    h1 = plot(1:length(errors{1}.grasp{jj}(:,1) ), errors{1}.grasp{jj}(:,1), 'rx');
    hold on
    h2 = plot(1:length(errors{1}.grasp{jj}(:,2) ), errors{1}.grasp{jj}(:,2), 'bo');
    hold on
end
hold off

strr = sprintf('Error values for tissue type = %s',enumTypes{1});
title(strr,'FontSize',fontS);
xlabel('time steps','FontSize',fontS)
ylabel('error','FontSize',fontS)
legend([h1,h2],'Correct Class','Incorrect Class')


figure(12)
for jj = 1:totalGrasps(2)
    
    h3 = plot(1:length(errors{2}.grasp{jj}(:,2) ), errors{2}.grasp{jj}(:,2), 'rx');
    hold on
    h4 = plot(1:length(errors{2}.grasp{jj}(:,1) ), errors{2}.grasp{jj}(:,1), 'bo');
    hold on
end

hold off

strr = sprintf('Error values for tissue type = %s',enumTypes{2});
title(strr,'FontSize',fontS);
xlabel('time steps','FontSize',fontS)
ylabel('error','FontSize',fontS)
legend([h3,h4],'Correct Class','Incorrect Class')


figure(13)
for jj = 1:totalGrasps(1)
    h1 = plot(1:length(error_sums{1}.grasp{jj}(:,1) ), error_sums{1}.grasp{jj}(:,1), 'r-');
    hold on
    h2 = plot(1:length(error_sums{1}.grasp{jj}(:,2) ), error_sums{1}.grasp{jj}(:,2), 'b--');
    hold on
end
hold off

strr = sprintf('Error sums over time for tissue type = %s',enumTypes{1});
title(strr,'FontSize',fontS);
xlabel('time steps','FontSize',fontS)
ylabel('error','FontSize',fontS)
legend([h1,h2],'Correct Class','Incorrect Class')


figure(14)
for jj = 1:totalGrasps(2)
    h3 = plot(1:length(error_sums{2}.grasp{jj}(:,2) ), error_sums{2}.grasp{jj}(:,2), 'r-');
    hold on
    h4 = plot(1:length(error_sums{2}.grasp{jj}(:,1) ), error_sums{2}.grasp{jj}(:,1), 'b--');
    hold on
end

hold off

strr = sprintf('Error sums over times for tissue type = %s',enumTypes{2});
title(strr,'FontSize',fontS);
xlabel('time steps','FontSize',fontS)
ylabel('error','FontSize',fontS)
legend([h3,h4],'Correct Class','Incorrect Class')

offset1 = 100000;
offset2 = 300000;

figure(15)
for jj = 1:totalGrasps(1)
    alphas = ( error_sums{1}.grasp{jj}(:,2) - error_sums{1}.grasp{jj}(:,1) ) ./ ( error_sums{1}.grasp{jj}(:,2) + offset1 );
    
    h5 = plot(1:length(alphas),alphas,'r-');
    hold on
end

hold off

strr = sprintf('Confidence over time for tissue type = %s',enumTypes{1});
title(strr,'FontSize',fontS);
xlabel('time steps','FontSize',fontS)
ylabel('alpha','FontSize',fontS)


figure(16)
for jj = 1:totalGrasps(2)
    alphas = ( error_sums{2}.grasp{jj}(:,1) - error_sums{2}.grasp{jj}(:,2) ) ./ (error_sums{2}.grasp{jj}(:,1) + offset2 );
    
    h6 = plot(1:length(alphas),alphas,'r-');
    hold on
end

hold off

strr = sprintf('Confidence over time for tissue type = %s',enumTypes{2});
title(strr,'FontSize',fontS);
xlabel('time steps','FontSize',fontS)
ylabel('alpha','FontSize',fontS)


