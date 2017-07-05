function h = windowDiscriminantPlotGrid(Data, Y, pthresh)

    if(nargin == 2)
       pthresh = 0.3; 
    end
    %Nice Grid For Plotting
    DIV = [100, 100];
    [NN,SS] = size(Data);
    classes = unique(Y);

    %making grid for doing stuff
    limz = [min(Data);max(Data)];
    Grid = ndimgrid(limz,DIV);

%     kernelsz = 1.06*norm(cov(Data))*NN^(-1/5) %bandwidth estimate
    kernelsz = mean(range(Data))*0.1; %range percentage;
%     kernelsz = min([bwest,rks])
%     kernelsz = mean(std(Data))

    for kk = 1:length(classes)
        datac{kk} = Data(Y == classes(kk),:);
    end

    minpoints = 10;

    for ii = 1:length(Grid)
        testpoint = Grid(ii,:);
        
        
        for kk = 1:length(classes)
            windowData = DataInWindow(datac{kk},testpoint,kernelsz);
           
            
            P(kk) = probabilityWindow(windowData,testpoint);
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
    probgrid(:,3) = probgrid(:,3) ./ max(probgrid(:,3));
    K = {};
    
    h = figure;
    for kk = 1:length(classes)
        scatter(datac{kk}(:,1),datac{kk}(:,2));
        hold on
        str = sprintf('class %d',kk);
        K{kk} = str;
    end
    handle = Surface3D(Grid(:,1),Grid(:,2),probgrid(:,3));
    hold off
    title('Seperability Windows','fontsize', 12)
    xlabel('x1','FontSize', 12)
    ylabel('x2','FontSize', 12)
    zlabel('sep','FontSize', 12)
    hl = legend(K,'Location','northeast');
    set(hl,'FontSize',12);
    az = -33;
    el = 20;
    view(az, el);
    
    
    %Only keep data if seperable wont work for more than binary
    subdata{1} = probgrid(probgrid(:,3) > pthresh ,[1 2]); %class 1
    subdata{2} = probgrid(probgrid(:,3) < -pthresh ,[1 2]); %class 2
    
    figure
    scatter(subdata{1}(:,1),subdata{1}(:,2),'r.')
    hold on
    scatter(subdata{2}(:,1),subdata{2}(:,2),'b.')
    hold on
    handle = Surface3D(Grid(:,1),Grid(:,2),probgrid(:,3));
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