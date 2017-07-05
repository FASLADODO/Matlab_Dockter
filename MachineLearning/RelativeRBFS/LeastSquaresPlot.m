function [] = LeastSquaresPlot(Data,Labels,LS_Param,modelfunc,cols)
%Data = class data to train on (rows are samples, columns are dimensions)
%Labels = column of unique class lables with same length as Data [1;2...]
%LS_Param: come from LeastSquaresClassification.m
% modelfunc eg: = @(X) [X(:,1),X(:,2),X(:,1).^2]
% outputcolumn: column for Y = Data(:,outputcolumn)

    [NN,SS] = size(Data);
    cslist = unique(Labels);
    
    linD = linspaceND(Data,NN);
    
    figure
    if(SS == 2)
        gscatter(Data(:,cols(1)),Data(:,cols(2)),Labels,'br');
    elseif(SS ==3)
        gscatter3(Data(:,cols(1)),Data(:,cols(2)),Data(:,cols(3)),Labels,'br');
    end
    
    %find least squares params
    for cc = 1:length(LS_Param)
        Dtemp = linD;
        %Get training info
        DOn = modelfunc(Dtemp);

        %calulcate least squares params
        Y_Estimate = DOn*LS_Param{cc}.Params;
        
        hold on
        if(SS == 2)
            scatter(linD(:,cols(1)),Y_Estimate,'g.');
        elseif(SS == 3)
            scatter3(linD(:,cols(1)),linD(:,cols(2)),Y_Estimate,'g.');
        end
    end
    hold off
    title('data with least squares fit')
    xlabel('x1')
    ylabel('x2')

end