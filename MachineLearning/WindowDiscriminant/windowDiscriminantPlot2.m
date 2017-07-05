function [h,Mdl] = windowDiscriminantPlot2(Data, Y, pthresh)

    %window data using K-nearest neighbors instead of window size

    kn = floor(0.01*length(Data))
    
    rod_is_dumb=false;

    if(nargin == 2)
       pthresh = 0.3; 
    end

    [NN,SS] = size(Data);
    classes = unique(Y);


    for kk = 1:length(classes)
        datac{kk} = Data(Y == classes(kk),:);
        Mdl{kk} = KDTreeSearcher(datac{kk});
        %Mdl{kk} = createns(datac{kk},'nsmethod','kdtree');
    end

    minpoints = 10;

    for ii = 1:length(Data)
        testpoint = Data(ii,:);
        
        if(rod_is_dumb)
            figure(1)
            scatter(testpoint(:,1),testpoint(:,2),'g*')
            hold on
        end
        
        for kk = 1:length(classes)
            %windowData = DataInWindow(datac{kk},testpoint,kernelsz);
            idwd = knnsearch(Mdl{kk},testpoint,'k',kn);
            windowData = Mdl{kk}.X(idwd,:);
            
            if(rod_is_dumb)
                if(kk == 1)
                    scatter(windowData(:,1),windowData(:,2),'r.')
                else
                    scatter(windowData(:,1),windowData(:,2),'b.')
                end
                hold on
            end
            
            P(kk) = probabilityWindow(windowData,testpoint);
        end
        if(rod_is_dumb)
            hold off
            pause(0.03)
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
    handle = Surface3D(probgrid(:,1),probgrid(:,2),probgrid(:,3));
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
    handle = Surface3D(probgrid(:,1),probgrid(:,2),probgrid(:,3));
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