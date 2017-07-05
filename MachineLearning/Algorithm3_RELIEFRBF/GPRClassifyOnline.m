function [Class,ClassTime,ScoreTime,muest] = GPRClassifyOnline(Model,XTest)
%online classification with GPR
% Model trained with GPRClassifyTrain.m
% data to be classified
% Returns:
% Class is the overall class estimate for the segment
% ClassTime is the per timestep class estimate
% ScoreTime is the running sum score at each time step

[NN,SS] = size(XTest);

%see if we should scale
if(~isempty(Model.shift) && ~isempty(Model.scale) )
    [XTest] = MeanVarianceScale(XTest,Model.shift,Model.scale);
end

%now lets get our 2 covariances (K_*N,K_**)
K_test_train = GaussianKernel(XTest,Model.XTrain,Model.kparams);
K_test_test = GaussianKernel(XTest,XTest,Model.kparams);

%we compute this only once since we use it twice
LK = K_test_train * Model.K_inv;

%figure out per time step stuff
muest = LK * Model.YTrain;
ScoreTime = cumsum(muest);
ClassTime = Model.f_classify(ScoreTime);
Class = Model.f_classify(ScoreTime(end));

%posterior variance (not used)
sigmaest = (K_test_test + Model.sigman*eye(NN) ) - LK* K_test_train';


end