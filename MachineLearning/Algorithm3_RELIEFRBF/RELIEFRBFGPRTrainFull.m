function [Model] = RELIEFRBFGPRTrainFull(Data,Labels,rw)
%performs full subsampling based on reliefrbf and train gpr classification
%model

%scale it
[Data,shift,scale] = MeanVarianceScale(Data);

%subsample the data
[subdiff,subdata,sublabels] = GPRSubsample(Data,Labels, rw);

%train gpr classification
Model = GPRClassifyTrain(subdata,sublabels);

%add in scaling
Model.shift = shift;
Model.scale = scale;

end