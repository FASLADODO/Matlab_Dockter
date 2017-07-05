function [Models, kmg] = windowDiscriminantTrain3D(Data, Y, pthresh, kmg, kmeanlimit, option)
    %Data is column vectors of data in each state
    %Y is column vector of classification labels
    % pthresh is the threshold below which data is considered inseperable
    %kmg is the number of gaussian mixture models (starts as undefined
    
    
    if(nargin == 4)
       option = 'plotoff'; 
    end

    [NN,SS] = size(Data);
    classes = unique(Y);

%     kernelsz = 1.06*norm(std(Data))*NN^(-1/5) %bandwidth estimate
    kernelsz =  mean(range(Data))*0.1; %range percentage;
%     kernelsz = min([bwest,rks])
%     kernelsz = min(sqrt(std(Data)))
%     kernelsz = mean(std(Data))


    for kk = 1:length(classes)
        %Get seperate data for each class
        datac{kk} = Data(Y == classes(kk),:);
        idxc{kk} = find(Y == classes(kk));
    end

    minpoints = 10;

    %loop through all training data points
    for ii = 1:length(Data)
        testpoint = Data(ii,:);

        %now get all data in a window for each class
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

        probwindowd(ii,:) = [ testpoint, onclass - offclass] ;
    end

    %Scale this
    pscale = 1 / max(abs(probwindowd(:,end)));
    probwindowd(:,end) = probwindowd(:,end) .*pscale;

    %Only keep data if seperable wont work for more than binary
    subdata{1} = probwindowd(probwindowd(:,end) > pthresh ,[1,2,3]); %class 1
    subdata{2} = probwindowd(probwindowd(:,end) < -pthresh ,[1,2,3]); %class 2
    
    subP{1} = probwindowd(probwindowd(:,end) > pthresh ,end); %class 1
    subP{2} = probwindowd(probwindowd(:,end) < -pthresh ,end); %class 2
    
    %Determine number of gaussian mixture components ?
    if(isempty(kmg))
        for cc = 1:length(classes)
            eva = evalclusters(subdata{cc},'kmeans','gap','KList',[1:kmeanlimit]);
            km(cc) = eva.OptimalK; %assume this is solved
        end
        kmg = km;
    else
       km = kmg; 
    end

    Models = [];
    for cc = 1:length(classes)
       dataon = subdata{cc};
       %Dope clustering function to find clumps of data in each class
       %Tclust = clusterdata(dataon,'linkage','centroid','maxclust',km(cc));
       [meansk,Tclust] = kMeansIterative(dataon,km(cc)); %Rods kmeans version
       for mm = 1:km(cc)
           idxtemp = find(Tclust == mm);
           Models{cc}.cluster{mm}.mu = meansk(mm,:); %mean(dataon(idxtemp,:));
           Models{cc}.cluster{mm}.sigma = cov(dataon(idxtemp,:));
           Models{cc}.cluster{mm}.scale = 1/gaussianProbMV(Models{cc}.cluster{mm}.mu,Models{cc}.cluster{mm}.sigma,Models{cc}.cluster{mm}.mu);
       end
       meansk
    end

    if(option == 'ploton')
        figure
        for kk = 1:length(classes)
            scatter3(subdata{kk}(:,1),subdata{kk}(:,2),subdata{kk}(:,3),20,subP{kk},'.');
            hold on
            str = sprintf('class %d',kk);
            K{kk} = str;
        end
        hold on
        cc = 1;
        for mm = 1:km(cc) 
            scatter3(Models{cc}.cluster{mm}.mu(1),Models{cc}.cluster{mm}.mu(2),Models{cc}.cluster{mm}.mu(3),100,'g*')
            hold on
        end
        cc = 2;
        for mm = 1:km(cc) 
            scatter3(Models{cc}.cluster{mm}.mu(1),Models{cc}.cluster{mm}.mu(2),Models{cc}.cluster{mm}.mu(3),100,'k*')
            hold on
        end

        hold off
        title('Seperated Data','FontSize', 12)
        xlabel('x1','FontSize', 12)
        ylabel('x2','FontSize', 12)
        hl = legend('class1','class2','Location','northeast');
        set(hl,'FontSize',12);
        colormap(cool);
        colorbar;
        az = -33;
        el = 20;
        view(az, el);
    end

end