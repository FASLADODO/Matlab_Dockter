%%
close all


%% importing raw data
clear all

endData = 21000;
countsPerRev = 500; %encoder (astrini thesis)

filename = {'bladder4Hz1.lvm';'bladder4Hz2.lvm';'gallbladder4Hz1.lvm';'gallbladder4Hz2.lvm';'liver4Hz1.lvm';'liver4Hz2.lvm'; 'smallbowel4Hz1.lvm';'smallbowel4Hz2.lvm'};
enumTypes = {'bladder'; 'gallbladder'; 'liver' ;'smallbowel' };
tissuetypeIn = [1,1,2,2,3,3,4,4];

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
    
    ax1 = subplot(3,1,1);
    plot(FullData{kk}.angle(1:plotlength));
    ylabel('angle')
    
    ax2 = subplot(3,1,2);
    plot(FullData{kk}.angledot(1:plotlength));
    ylabel('angledot')
    
%     ax3 = subplot(4,1,3)
%     plot(FullData{kk}.angledotdot(1:plotlength))
%     ylabel('angledotdot')
%     
    ax3 = subplot(3,1,3);
    plot(FullData{kk}.stress(1:plotlength));
    ylabel('stress')

    
    linkaxes([ax1,ax2,ax3],'x')
    
    str = sprintf('Raw Data: Grasp #: %i, tissue: %s',kk, enumTypes{ FullData{kk}.tissuetype } );
    suptitle(str);
    
end

%% GO BACKWARD IN TIME!

defgraspstep = 0.15;

% SIGNS ARE CORRECT FOR RAW
stressThreshold = 1.3; %Newtons, Critical threshold at which to start counting a grasp. As small as possible!
stressThreshHigh = 3; %Newtons, Critical threshold at which to start counting a grasp. As small as possible!
VelThreshold = 1; %mm/s, minimum positive closing speed at which to consider grasp as occuring
MinGraspInterval = 30; %Minimum number of time steps between grasps
MinGraspStep = defgraspstep; %Minimum namount of time a grasp should last
MinAngleDisp = 0.9;

velslow = 0.5;

sumlen = 0;
inder = 0;

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
                if (FullData{kk}.t(ss) < FullData{kk}.t(GraspEnd_i) - MinGraspStep) && (FullData{kk}.stress(ss) < stressThreshold) && (FullData{kk}.angle(GraspEnd_i) - FullData{kk}.angle(ss) > MinAngleDisp && FullData{kk}.angledot(ss) - FullData{kk}.angledot(ss+1) < velslow)
                    
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
                    
                    
                    %sumlen = sumlen + length(segDataTemp{kk}.GraspData{nGrasp}.stress);
                    sumlen = sumlen + (segDataTemp{kk}.GraspData{nGrasp}.t(end) - segDataTemp{kk}.GraspData{nGrasp}.t(1));
                    inder = inder + 1;
                    
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
    segData{kk}.strain = [];
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
        segData{kk}.strain = [segData{kk}.strain; (segData{kk}.GraspData{indnew}.angle - segData{kk}.GraspData{indnew}.angle(1))./segData{kk}.GraspData{indnew}.angle(1) ];
        
        fprintf('1: %f, end: %f \n',segData{kk}.GraspData{indnew}.angle(1),segData{kk}.GraspData{indnew}.angle(end));
        
        segData{kk}.angle = [segData{kk}.angle; segData{kk}.GraspData{indnew}.angle];
        segData{kk}.angledot = [segData{kk}.angledot; segData{kk}.GraspData{indnew}.angledot];
        segData{kk}.angledotdot = [segData{kk}.angledotdot; segData{kk}.GraspData{indnew}.angledotdot];
        segData{kk}.tissuetype = FullData{kk}.tissuetype;
        
        indnew = indnew + 1;
    end
    
    figure(kk)
    ax1 = subplot(3,1,1);
    plot(FullData{kk}.t,FullData{kk}.angle);
    hold on
    scatter(segData{kk}.t,segData{kk}.angle,'.r');
    hold off
    ylabel('angle (rad)')
    
    ax2 = subplot(3,1,2);
    plot(FullData{kk}.t,FullData{kk}.angledot);
    hold on
    scatter(segData{kk}.t,segData{kk}.angledot,'.r');
    hold off
    ylabel('angledot (rad)')
    
    ax3 = subplot(3,1,3);
    plot(FullData{kk}.t,FullData{kk}.stress);
    hold on
    scatter(segData{kk}.t,segData{kk}.stress,'.r');
    hold off
    ylabel('stress (N)')
    
    linkaxes([ax1,ax2,ax3],'x')
    
    str = sprintf('Segmented Data: Grasp #: %i, tissue: %s',kk, enumTypes{ FullData{kk}.tissuetype } );
    suptitle(str);
    
end

sumlen/inder


disp('Number of grasps:')
for ii = 1:max(tissuetypeIn)
    str = sprintf('tissue type: %s, # grasps: %i \n',enumTypes{ii},length(segData{ii}.GraspData));
    disp(str)
end


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
CM = { [1,0,0] , [0,1,0], [0,0,1], [0,1,1], [0,1,1], [1,0,1], [0,0,0], [.5,1,1] };


figure

for kk = 1:numFiles
    h(tissuetypeIn(kk)) = scatter3(segData{kk}.angle,segData{kk}.angledot,segData{kk}.stress,'.','MarkerEdgeColor',CM{FullData{kk}.tissuetype});
    
    %legend_name{kk} = enumTypes{ FullData{kk}.tissuetype };
    hold on
end
hold off

title('Phase Portrait')
xlabel('angle (rad)')
ylabel('Angledot (rad/s)')
zlabel('strss (N)')
legend([h(1),h(2),h(3),h(4)],enumTypes{:})


figure

for kk = 1:numFiles
    h(tissuetypeIn(kk)) = scatter(segData{kk}.strain,segData{kk}.angledot,'.','MarkerEdgeColor',CM{FullData{kk}.tissuetype});
    
    %legend_name{kk} = enumTypes{ FullData{kk}.tissuetype };
    hold on
end
hold off

title('Phase Portrait (strain, angledot)')
xlabel('strain (rad)')
ylabel('Angledot (rad/s)')
legend([h(1),h(2),h(3),h(4)],enumTypes{:})


figure
for kk = 1:numFiles
    h(tissuetypeIn(kk)) = quiver(segData{kk}.angle,segData{kk}.angledot,segData{kk}.angledot,segData{kk}.angledotdot,'color',CM{FullData{kk}.tissuetype});
    
    %legend_name{kk} = enumTypes{ FullData{kk}.tissuetype };
    hold on
end
hold off

title('Phase Portrait (angle, angledot)')
xlabel('angle (rad)')
ylabel('Angledot (rad/s)')
legend([h(1),h(2),h(3),h(4)],enumTypes{:})


figure
for kk = 1:numFiles
    h(tissuetypeIn(kk)) = quiver(segData{kk}.angledot,segData{kk}.angledotdot,segData{kk}.angledotdot,segData{kk}.angledotdot,'color',CM{FullData{kk}.tissuetype});
    
    %legend_name{kk} = enumTypes{ FullData{kk}.tissuetype };
    hold on
end
hold off

title('Phase Portrait (angledot, angledotdot)')
xlabel('Angledot (rad/s)')
ylabel('Angledotdot (rad/s^2)')
legend([h(1),h(2),h(3),h(4)],enumTypes{:})

figure

for kk = 1:numFiles
    f(tissuetypeIn(kk)) = quiver(segData{kk}.strain,segData{kk}.stress,segData{kk}.angledot,segData{kk}.stressdot,'color',CM{FullData{kk}.tissuetype});
    
    %legend_name{kk} = enumTypes{ FullData{kk}.tissuetype };
    hold on
end
hold off

title('Phase Portrait (strain, stress)')
xlabel('strain ()')
ylabel('stress (N)')
legend([f(1),f(2),f(3),f(4)],enumTypes{:})


figure

for kk = 1:numFiles
    f(tissuetypeIn(kk)) = quiver(segData{kk}.stress,segData{kk}.angledot,segData{kk}.stressdot,segData{kk}.angledotdot,'color',CM{FullData{kk}.tissuetype});
    
    %legend_name{kk} = enumTypes{ FullData{kk}.tissuetype };
    hold on
end
hold off

title('Phase Portrait files (stress, angledot)')
xlabel('stress (N)')
ylabel('angledot (rad/s)')
legend([f(1),f(2),f(3),f(4)],enumTypes{:})



%% combine segments

%set to empty
for jj = 1:max(tissuetypeIn)
    combineSeg{jj}.t = [];
    combineSeg{jj}.stress = [];
    combineSeg{jj}.stressdot = [];
    combineSeg{jj}.stressdotdot = [];
    combineSeg{jj}.strain = [];
    combineSeg{jj}.angle = [];
    combineSeg{jj}.angledot = [];
    combineSeg{jj}.angledotdot = [];
    combineSeg{jj}.D = [];
    combineSeg{jj}.U = [];
    combineSeg{jj}.Phi = [];
end

% combined segmented grasps
for kk = 1:numFiles
    combineSeg{FullData{kk}.tissuetype}.t = [combineSeg{FullData{kk}.tissuetype}.t ; segData{kk}.t ];
    combineSeg{FullData{kk}.tissuetype}.stress = [combineSeg{FullData{kk}.tissuetype}.stress ; segData{kk}.stress ];
    combineSeg{FullData{kk}.tissuetype}.stressdot = [combineSeg{FullData{kk}.tissuetype}.stressdot ; segData{kk}.stressdot ];
    combineSeg{FullData{kk}.tissuetype}.stressdotdot = [combineSeg{FullData{kk}.tissuetype}.stressdotdot ; segData{kk}.stressdotdot ];
    combineSeg{FullData{kk}.tissuetype}.strain = [combineSeg{FullData{kk}.tissuetype}.strain ; segData{kk}.strain ];
    
    combineSeg{FullData{kk}.tissuetype}.angle = [combineSeg{FullData{kk}.tissuetype}.angle ; segData{kk}.angle ];
    combineSeg{FullData{kk}.tissuetype}.angledot = [combineSeg{FullData{kk}.tissuetype}.angledot ; segData{kk}.angledot ];
    combineSeg{FullData{kk}.tissuetype}.angledotdot = [combineSeg{FullData{kk}.tissuetype}.angledotdot ; segData{kk}.angledotdot ];
    
end

%Create D and U matrix for RLS
for jj = 1:max(tissuetypeIn)
    %combineSeg{jj}.D = [combineSeg{jj}.angledotdot,combineSeg{jj}.angledotdot.^(1/2),combineSeg{jj}.angledotdot.^2,combineSeg{jj}.angledotdot.^3, combineSeg{jj}.angledot, combineSeg{jj}.angledot.^(1/2),combineSeg{jj}.angledot.^2, combineSeg{jj}.angledot.^3, combineSeg{jj}.angle, combineSeg{jj}.angle.^(1/2),combineSeg{jj}.angle.^2, combineSeg{jj}.angle.^3, ones(length(combineSeg{jj}.angle) , 1) ];
    combineSeg{jj}.D = [combineSeg{jj}.angledotdot,combineSeg{jj}.angledotdot.^2,combineSeg{jj}.angledotdot.^3, combineSeg{jj}.angledot,combineSeg{jj}.angledot.^2, combineSeg{jj}.angledot.^3, combineSeg{jj}.angle,combineSeg{jj}.angle.^2, combineSeg{jj}.angle.^3, ones(length(combineSeg{jj}.angle) , 1) ];
    
    combineSeg{jj}.U = [combineSeg{jj}.stress ];
    combineSeg{jj}.Phi = inv(combineSeg{jj}.D'*combineSeg{jj}.D) * combineSeg{jj}.D' * combineSeg{jj}.U;
end


%Akaike information criteria
[NN,SS] = size(combineSeg{1}.D);

states = {'All','no dynamics'}; %,'only dynamics'};
% labels = {'xdotdot','xdotdot^(1/2)','xdotdot^2','xdotdot^3','xdot','xdot^(1/2)','xdot^2','xdot^3','x','x^(1/2)','x^2','x^3','1'}; %this affects model in
labels = {'thetaddot','thetaddot^2','thetaddot^3','thetadot','thetadot^2','thetadot^3','theta','theta^2','theta^3','1'}; %this affects model in

for jj = 1:max(tissuetypeIn)
        X_in = combineSeg{jj}.D;
        Y_in = combineSeg{jj}.U;
        [AIC{jj},bestcombos{jj}] = AkaikeInformationCriteria(Y_in,X_in);
        
        minner = [];
        %plot that shit
        figure
        for kk = 1:SS
            minner(kk) = min(nonzeros(AIC{1}(:,kk)));
            dd = nonzeros(AIC{jj}(:,kk));
            h = scatter(ones(length(dd),1)*kk,dd);
            hold on
        end
        hold off
        %set(gca,'XTickLabel',states)
        xlabel('# Features')
        ylabel('AIC')
        sttit = sprintf('Akaike Information: %s',enumTypes{jj});
        title(sttit)
        
        figure
        diffin = abs(diff(minner));
        plot(1:(SS-1),diffin);
        xlabel('# Features')
        ylabel('Min AIC')
        sttit = sprintf('AIC Min Diff: %s',enumTypes{jj});
        title(sttit)
        
        %print off best at 7
        sb = 8;
        arr = bestcombos{jj}{sb};
        besties{jj} = labels(arr);
end


% %DO SOME BAYESIAN INFORMATION STUFF!
% states = {'xdotdot','xdot','x','x^2','x^3','1'};
% for jj = 1:max(tissuetypeIn)
%     X_in = combineSeg{jj}.D;
%     Y_in = combineSeg{jj}.U;
%     [Errorall{jj}, MeanErrorall{jj}, Paramsall{jj}] = BayesianInformationCriteria(X_in,Y_in,states,enumTypes{jj},'TLS');
% end


%should evaluate to:
phi1 = [-0.000621406370117322;0.0976101250736223;-0.388256674969638;0.131870191736068;0.0373350034929734];
phi2 = [-0.000826176494778239;0.101152437533438;0.503038080455392;-0.429175896084507;0.129236022201036];
phi3 = [-0.00174977066715285;0.131433006396432;-0.246563945926546;0.0476950832392320;0.0553277583906477];
phi4 = [-0.000617185498974871;0.0233326316387804;0.679253921276928;-0.238591813349058;0.0703108547058125];


%% Online leave one out

ind = 1;
classify = [];
%loop through all files

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
                        online.state = [segData{ii}.GraspData{hh}.angledotdot, segData{ii}.GraspData{hh}.angledot, segData{ii}.GraspData{hh}.angle, segData{ii}.GraspData{hh}.angle.^2, segData{ii}.GraspData{hh}.angle.^3, ones(length(segData{ii}.GraspData{hh}.angle) ,1 ) ];
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
        
        %get parameters for all tissue types
        for qq = 1:max(tissuetypeIn)
            trainSeg{qq}.Phi = inv(trainSeg{qq}.state'*trainSeg{qq}.state) * trainSeg{qq}.state' * trainSeg{qq}.input;
        end
        
        %test parameters for all tissue types
        error_temp = [];
        for qq = 1:max(tissuetypeIn)
            anss = cumsum(abs(online.input - online.state*trainSeg{qq}.Phi ));
            error_temp = [error_temp, anss];
        end
        
        %Get the minimum sum error
        [minner,idx] = min(error_temp(end,:));
        %save actual class and estimated class
        if true % (segData{kk}.tissuetype ~= 4)
            classify(ind,:) = [FullData{kk}.tissuetype, idx];
            ind = ind + 1;
        end
    end
end

%Get classification precentage
correct = (classify(:,1) == classify(:,2) );
disp('classification percentage: ')
corr_percent = sum(correct)/length(correct)