function [model,Diff] = RelativeGaussianTrain(Data,Input,processfunc,relativefunc,option,beta0)

%fit parameters to a linear or nonlinear model
%then train a KL divergence fit to get relative seperability
%Data: data matrix, rows as samples
%Input: column vector eg forces
% if option='linear'
% processfunc = @(X) ([X(:,1),X(:,2),X(:3,)]) columns for process fit
% relativefunc = @(X) ([X(:,1),X(:,2)]) columns for process fit
% if option='nonlinear'
% processfunc = = @(Param,X) Param(1)*X(:,1) + Param(2)*X(:,2) + Param(3)*exp(-Param(4)*X(:,3) ) columns for nlinfit
% relativefunc = @(X) ([X(:,1),X(:,2)]) columns for process fit
%beta0 = param estimates [params1; params2]

%sizes
cslist = length(Data);
[NN,SS] = size(Data{1});

ParamTrue = [];
%Get fit params
for gg = 1:cslist
    dtemp = processfunc(Data{gg});
    utemp = Input{gg};
    
    if(strcmp(option,'linear'))
        ParamTrue{gg} = pinv(dtemp)*utemp;
    elseif(strcmp(option,'nonlinear'))
        ParamTrue{gg} = nlinfit(dtemp,utemp,processfunc,beta0(gg,:));
    end
end

Diff = [];
GaussianProbs = [];
%get gaussian from process
for gg = 1:cslist
    dtemp = processfunc(Data{gg});
    utemp = Input{gg};
    
    %test classification with all models
    for cc = 1:cslist
        if(strcmp(option,'linear'))
            [PL,~] = LinearGaussianSimple(dtemp,utemp,ParamTrue{cc});
        elseif(strcmp(option,'nonlinear'))
            [PL,~] = NonLinearGaussian(dtemp,utemp,ParamTrue{cc},processfunc);
        end
        
        if(gg == cc)
            pwithin = PL;
        else
            pbetween = PL;
        end
    end
    Diff = [Diff; ComputeRBFDifference(pwithin,pbetween)];
    GaussianProbs = [GaussianProbs; pwithin];
end

% Now train the relative seperability to some new params
DataAll = [];
for gg = 1:cslist
    dtemp = Data{gg};
    DataAll = [DataAll; dtemp];
end
%train and get scale
drel = relativefunc(DataAll);
param_diff = pinv(drel)*Diff;
param_diff_scale = max(Diff);
   

%stash in struct
model.ProcessParams = ParamTrue;
model.SeperableParams = param_diff;
model.SeperableParamScale = param_diff_scale;
model.ProcessFunction = processfunc;
model.SeperableFunction = processfunc;
model.option = option;  
    
end


