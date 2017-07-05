function h = windowDiscriminantPlot3D(Data, Y, pthresh)

    if(nargin == 2)
       pthresh = 0.3; 
    end
    %Nice Grid For Plotting
    DIV = [100, 100];
    [NN,SS] = size(Data);
    classes = unique(Y);

%     kernelsz = 1.06*norm(cov(Data))*NN^(-1/5) %bandwidth estimate
    kernelsz = min(range(Data))*0.1 %range percentage;
%     kernelsz = min([bwest,rks])
%     kernelsz = mean(std(Data))

    for kk = 1:length(classes)
        datac{kk} = Data(Y == classes(kk),:);
        idxc{kk} = find(Y == classes(kk));
    end

    minpoints = 10;

    for ii = 1:length(Data)
        testpoint = Data(ii,:);
        
        
        for kk = 1:length(classes)
            windowData = DataInWindow(datac{kk},testpoint,kernelsz);
           
            
            if(length(windowData) < minpoints)
               P(kk) = 0; 
            else
                P(kk) = probabilityWindow2(windowData,testpoint);
            end
        end
        

        [maxr,idxm] = max(P);
        onclass = P(1);
        offclass = P(2);
        %offclass(idxm) = [];

        probgrid(ii,:) = [testpoint, onclass - offclass ];
    end

    %Scale this
    probgrid(:,end) = probgrid(:,end) ./ max(probgrid(:,end));
    K = {};
    
    h = figure;
    for kk = 1:length(classes)
        scatter3(datac{kk}(:,1),datac{kk}(:,2),datac{kk}(:,3),20,probgrid(idxc{kk},end),'.');
        hold on
        str = sprintf('class %d',kk);
        K{kk} = str;
    end

    hold off
    title('Seperability Windows','fontsize', 12)
    xlabel('x1','FontSize', 12)
    ylabel('x2','FontSize', 12)
    zlabel('sep','FontSize', 12)
    hl = legend(K,'Location','northeast');
    set(hl,'FontSize',12);
    colormap(cool);
    colorbar;
    az = -33;
    el = 20;
    view(az, el);
    
    
    %Only keep data if seperable wont work for more than binary
    subdata{1} = probgrid(probgrid(:,end) > pthresh ,[1 2 3]); %class 1
    subdata{2} = probgrid(probgrid(:,end) < -pthresh ,[1 2 3]); %class 2
    
    figure
    scatter3(subdata{1}(:,1),subdata{1}(:,2),subdata{1}(:,3),'r.')
    hold on
    scatter3(subdata{2}(:,1),subdata{2}(:,2),subdata{2}(:,3),'b.')
    hold off
    title('Seperated Data','FontSize', 12)
    xlabel('x1','FontSize', 12)
    ylabel('x2','FontSize', 12)
    hl = legend('class1','class2','Location','northeast');
    set(hl,'FontSize',12);
    az = -33;
    el = 20;
    view(az, el);


end