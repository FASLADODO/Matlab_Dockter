function [h,probgrid] = windowDiscriminantAlt(Data, Y)
    %Nice Grid For Plotting
    DIV = [30, 30];
    [NN,SS] = size(Data);
    classes = unique(Y);

    %making grid for doing stuff
    limz = [min(Data);max(Data)];
    Grid = ndimgrid(limz,DIV);

%     bwest = 1.06*norm(cov(Data))*NN^(-1/5); %bandwidth estimate
%     rks = range(Data)*0.1; %range percentage;
%     kernelsz = min([bwest,rks])
    kernelsz = min(sqrt(std(Data)))

    for kk = 1:length(classes)
        datac{kk} = Data(Y == classes(kk),:);
    end

    minpoints = 2*SS;

    for ii = 1:length(Data)
        testpoint = Data(ii,:);

        for kk = 1:length(classes)
            tempwd = DataInWindow(datac{kk},testpoint,kernelsz);
            if(length(tempwd) < minpoints)
               windowP{kk} = zeros(SS,1)*eps; 
            else
                windowP{kk} = probabilityWindowarray(tempwd);
            end
        end
        
        [ H_X, H_Y, H_XY ] = JointEntropy(windowP{1},windowP{2});

        probgrid(ii,:) = [testpoint, H_XY];
        
    end
    
    %Scale this
    probgrid(:,3) = probgrid(:,3) ./ max(probgrid(:,3));
    
    %Only keep data if seperable wont work for more than binary
    subdata = probgrid(probgrid(:,3) < 0.5 ,[1,2]); %class 1

    
    K = {};
    
    h = figure;
    scatter(subdata(:,1),subdata(:,2));
    hold on
    handle = Surface3D(Data(:,1),Data(:,2),probgrid(:,3));
    hold off
    title('Seperability Windows','fontsize', 12)
    xlabel('x1','FontSize', 12)
    ylabel('x2','FontSize', 12)
    zlabel('sep','FontSize', 12)
    az = -33;
    el = 20;
    view(az, el);


end