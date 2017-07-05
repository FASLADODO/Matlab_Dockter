function [] = RBFNicePlots(DataClass,Difference,Data_All,Class_Labels,columnlabels)

    [NN,SS] = size(Data_All);
    cs = length(DataClass);
    
    plotsize3 = 10;

    if(SS == 1)
        figure
        gscatter(Data_All,zeros(length(Data_All),1),Class_Labels,'br');
        for cc = 1:cs
            hold on
            scatter(DataClass{cc},Difference{cc},'g.');
        end
        hold off
        title('Seperability w/ Discriminant')
        xlabel(columnlabels(1))
        ylabel('Relative P')
        
        figure
        gscatter(Data_All,zeros(length(Data_All),1),Class_Labels,'br');
        xlabel(columnlabels(1))
    elseif(SS == 2)
        figure
        gscatter(Data_All(:,1),Data_All(:,2),Class_Labels,'br');
        for cc = 1:cs
            hold on
            Surface3D(DataClass{cc}(:,1),DataClass{cc}(:,2),Difference{cc},'mesh');
        end
        hold off
        title('Seperability w/ Discriminant')
        xlabel(columnlabels(1))
        ylabel(columnlabels(2))
        zlabel('Seperability')
        view(45, 45);
    elseif(SS == 3)
        figure
        for cc = 1:cs
            scatter3(DataClass{cc}(:,1),DataClass{cc}(:,2),DataClass{cc}(:,3),plotsize3,Difference{cc});
            hold on
        end
        title('Seperability w/ Discriminant')
        xlabel(columnlabels(1))
        ylabel(columnlabels(2))
        zlabel(columnlabels(3))
        colormap cool
        colorbar
        view(45, 45);
        
        figure
        gscatter3(Data_All(:,1),Data_All(:,2),Data_All(:,3),Class_Labels,'br');
        xlabel(columnlabels(1))
        ylabel(columnlabels(2))
        zlabel(columnlabels(3))
        view(45, 45);
        
    else
        warning('cannot plot this many dimensions');
        for jj = 2:SS
            figure
            gscatter(Data_All(:,jj-1),Data_All(:,jj),Class_Labels,'br');
            hold on
            for cc = 1:cs
                Surface3D(DataClass{cc}(:,jj-1),DataClass{cc}(:,jj),Difference{cc},'mesh');
                hold on
            end
            hold off
            title('Seperability w/ Discriminant')
            xlabel(columnlabels(jj-1))
            ylabel(columnlabels(jj))
            zlabel('sep')
            colormap cool
            colorbar
            view(45, 45);
        end
        
        
    end

end