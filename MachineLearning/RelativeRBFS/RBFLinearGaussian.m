function [LGModel] = RBFLinearGaussian(Model,modelfunc,outputcolumn)
%Take cluster data and then train linear gaussians
    
    LGModel.function = modelfunc;
    figure
    %Now train linear gaussians with that data
    for cc = 1:length(Model)
        for idc = 1:Model{cc}.TotalClusters
            Dtemp = Model{cc}.cluster{idc}.clusterdata;

            %Get linear dimension columns
            DRaw = Dtemp;
            DRaw(:,outputcolumn) = [];
            DTrain = modelfunc(DRaw);

            %get output variable
            YTrain = Dtemp(:,outputcolumn);

            scatter(DRaw(:,1),YTrain)
            hold on

            %model
            ModelLSG = LinearGaussianTrain(DTrain,YTrain);
            LGModel.Class{cc}.cluster{idc}.Model = ModelLSG;

            [PL,Y_Est] = LinearGaussianOnline(DTrain,YTrain,ModelLSG);

            plot(DRaw(:,1),Y_Est,'g*')
            hold on
            handle2 = Surface3D(DRaw(:,1),YTrain,PL);
            hold on
        end
    end
    hold off
    title('linear gaussian fit')
end