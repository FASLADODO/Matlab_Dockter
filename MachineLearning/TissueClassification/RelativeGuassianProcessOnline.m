function [Difference,Difference_Scaled] = RelativeGuassianProcessOnline(Data,model)
%Train parameters and scale to map Data to Difference
%model: comes from RelativeGuassianProcessTrain()

    %figure out equation to map theta, thetadot to difference
    dtemp = model.Function(Data);
    Difference = dtemp*model.Params;
    Difference_Scaled = Difference./model.Scale;

end