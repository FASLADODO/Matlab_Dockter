function [ weights,region,x1grid, x2grid ] = DPP_Discriminant_Train_V2( tData, window, pStat, plotOn )
%Discriminant Phase Portrait Training
% tData is the training data for all classes and all runs
% kk = number of classes
% ii = number of states
% nn = number of data points in all training data
% tData{kk}.state{ii}() = [1 x nn]
% window is the data division size
% pStat is the state index to use (eg [1,2,3] )
% plotOn is 0 or 1 for plots on

%plot font
fontS = 14;

%number of class
numClass = length(tData);

% find limits
minx1 = min(min(tData{1}.state{pStat(1)}),min(tData{2}.state{pStat(1)}));
maxx1 = max(max(tData{1}.state{pStat(1)}),max(tData{2}.state{pStat(1)}));
minx2 = min(min(tData{1}.state{pStat(2)}),min(tData{2}.state{pStat(2)}));
maxx2 = max(max(tData{1}.state{pStat(2)}),max(tData{2}.state{pStat(2)}));

% get grid
tickx1 = (maxx1 - minx1)/window;
tickx2 = (maxx2 - minx2)/window;
x1grid = linspace(minx1,maxx1,window+1);
x2grid = linspace(minx1,maxx2,window+1);

if(plotOn)
    % Horizontal grid 
    for k = 1:length(x2grid)
      line([x1grid(1) x1grid(end)], [x2grid(k) x2grid(k)])
    end
    hold on
    % Vertical grid
    for k = 1:length(x2grid)
      line([x1grid(k) x1grid(k)], [x2grid(1) x2grid(end)])
    end
    hold on
end


%plot means
plotdex = {'bx';'kx';'cx';'ox'};
for kk = 1:numClass
    for ii = 1:length(x1grid)-1
       for jj = 1:length(x2grid)-1
           %segment region
          i1 = find( (tData{kk}.state{pStat(1)} >= x1grid(ii) & tData{kk}.state{pStat(1)} <= x1grid(ii+1)) & (tData{kk}.state{pStat(2)} >= x2grid(jj) & tData{kk}.state{pStat(2)} <= x2grid(jj+1)) );
          indexer{ii}.ind{jj} = i1;
          
          %find mean in each region
          region{kk}.hor{ii}.vert{jj}.state{pStat(1)} = tData{kk}.state{pStat(1)}(i1);
          region{kk}.hor{ii}.vert{jj}.avg_state{pStat(1)} = mean(tData{kk}.state{pStat(1)}(i1));
          region{kk}.hor{ii}.vert{jj}.state{pStat(2)} = tData{kk}.state{pStat(2)}(i1);
          region{kk}.hor{ii}.vert{jj}.avg_state{pStat(2)} = mean(tData{kk}.state{pStat(2)}(i1));
          
          if(plotOn)
              if(~isempty(region{kk}.hor{ii}.vert{jj}.avg_state{pStat(1)}) & ~isempty(region{kk}.hor{ii}.vert{jj}.avg_state{pStat(2)}) )
                plot(region{kk}.hor{ii}.vert{jj}.avg_state{pStat(1)},region{kk}.hor{ii}.vert{jj}.avg_state{pStat(2)},plotdex{kk},'markersize',10);
                hold on
              end
          end
       end
    end
end


%Get entropy weights
count = 1;
lt = 100; %length threshold
avgW = 1;
for ii = 1:length(x1grid)-1
    for jj = 1:length(x2grid)-1
        bc_true = false;
        wc_true = false;
        %Check that all classes have sufficient data
        checkAll = 0;
        for(kk = 1:numClass)
            if(length(region{kk}.hor{ii}.vert{jj}.state{pStat(1)}) > lt && length(region{kk}.hor{ii}.vert{jj}.state{pStat(2)}) > lt)
                checkAll = checkAll + 1;
            end
        end
        if(checkAll == numClass)
            %get between class
            data = [];
            for(kk = 1:numClass)
                data = [data; region{kk}.hor{ii}.vert{jj}.state{pStat(1)},region{kk}.hor{ii}.vert{jj}.state{pStat(2)}];
            end
            var_bc = var(data);
            bc_true = true;
            
            %region{1}.hor{ii}.vert{jj}.X_bc = X_bc;
            %region{1}.hor{ii}.vert{jj}.Y_bc = Y_bc;
            region{1}.hor{ii}.vert{jj}.var_bc = var_bc;
            
            %Get within class
            var_wc = 0;
            for(kk = 1:numClass)
                data_wc = [region{kk}.hor{ii}.vert{jj}.state{pStat(1)},region{kk}.hor{ii}.vert{jj}.state{pStat(2)}];
                var_wc = var_wc + var(data_wc);

                %region{kk}.hor{ii}.vert{jj}.X_wc = X_wc1;
                %region{kk}.hor{ii}.vert{jj}.Y_wc = Y_wc1;
                region{kk}.hor{ii}.vert{jj}.var_wc  = var_wc ;
            end
            wc_true = true;
            
            %Save weights
            if(bc_true && wc_true)
                weights(ii,jj) = sum(var_bc) / (sum(var_wc)/2);
                avgW = (avgW + abs(sum(var_wc) / sum(var_bc)) ) / (2);
            else
                weights(ii,jj) = 0;
            end
            count = count + 1;
        else
            for(kk = 1:numClass)
                if(length(region{kk}.hor{ii}.vert{jj}.state{pStat(1)}) > lt && length(region{kk}.hor{ii}.vert{jj}.state{pStat(2)}) > lt )
                    weights(ii,jj) = avgW;
                    break
                else
                    weights(ii,jj) = 0;
                end
                region{1}.hor{ii}.vert{jj}.X_bc = [];
                region{1}.hor{ii}.vert{jj}.Y_bc = [];
                region{1}.hor{ii}.vert{jj}.var_bc = [];
                region{1}.hor{ii}.vert{jj}.X_wc = [];
                region{1}.hor{ii}.vert{jj}.Y_wc = [];
                region{1}.hor{ii}.vert{jj}.var_wc = [];
            end
            
        end
        
    end
end

% Scale weights 0 -1
weights(weights>0) = weights(weights>0) - min(weights(weights>0));
weights = weights / ( max(max(weights)));

if(plotOn)
    cc = 1;
    for ii = 1:length(x1grid)-1
        for jj = 1:length(x2grid)-1
%             x1plot(cc) = x1grid(ii) + (tickx1/2);
%             x2plot(cc) = x2grid(jj) + (tickx2/2);
            if(~isnan(region{1}.hor{ii}.vert{jj}.avg_state{pStat(1)}) && ~isnan(region{2}.hor{ii}.vert{jj}.avg_state{pStat(1)}) )
                x1plot(cc) = (region{1}.hor{ii}.vert{jj}.avg_state{pStat(1)} + region{2}.hor{ii}.vert{jj}.avg_state{pStat(1)})/2;
                x2plot(cc) = (region{1}.hor{ii}.vert{jj}.avg_state{pStat(2)} + region{2}.hor{ii}.vert{jj}.avg_state{pStat(2)})/2;
            else
                x1plot(cc) = x1grid(ii) + (tickx1/2);
                x2plot(cc) = x2grid(jj) + (tickx2/2);
            end
            he(cc) = weights(ii,jj);
            cc= cc +1;
        end
    end

    x1lin = linspace(min(x1plot),max(x1plot),33);
    x2lin = linspace(min(x2plot),max(x2plot),33);

    [X,Y] = meshgrid(x1lin,x2lin);
    f = scatteredInterpolant(x1plot',x2plot',he');
    Z = f(X,Y);
    
     Z(Z<0)=0;

    mesher = mesh(X,Y,Z); %interpolated

    hold off
    view(0, 20);
    set(mesher,'facecolor','none')
    title('Phase Portrait with seperability weights','FontSize',fontS)
    xlabel('Angle (rad)','FontSize',fontS)
    ylabel('Angledot (rad/s)','FontSize',fontS)
    zlabel('Varaince Weights','FontSize',fontS)
end



end


