function Model = GPRClassifyTrain(X,Y,sigman,kparams)
%train a gaussian process classification model using data
%https://pdfs.semanticscholar.org/9b21/6ab14159c9193d6b10ddb370c636e9557c67.pdf
%https://www.quora.com/How-can-I-use-Gaussian-processes-to-perform-regression
%https://www.cs.toronto.edu/~hinton/csc2515/notes/gp_slides_fall08.pdf
%https://en.wikipedia.org/wiki/Gaussian_process

%X data matrix, rows as samples
%Y output column vector
%sigman: 0.5 observation noise
%kparams: kernel parameters = [lengthscale, sigmaf] sigmaf: %process noise
%lengthscale: typical distance between peaks
%A conditional of a gaussian distribution is also gaussian

%check args
if(nargin < 4)
    sigman = std(Y); %observation noise
    sigmaf = mean(std(X)); %process noise
    %silvermans rule of thumb
    lscale = 1.06 * sigman * length(Y)^(-1/5) ; %length scale
    kparams = [lscale, sigmaf];
end

Model = [];

%data size
[NN,SS] = size(X);

%grab all of the data to train with
XTrain = X; %used in model creation
YTrain = Y; %used in model creation

% now we compute our covariances, kernels, and means
%now lets get our 1 covariance (K_N)
K_train_train = GaussianKernel(XTrain,XTrain,kparams);

Kns = K_train_train + sigman*eye(NN);

if(false)
    %We're gonna optimize our length scale
    %marginal likelihood (how good are our current params)
    plog_prev = -0.5*YTrain'*(Kns\YTrain) - 0.5*log10(det(Kns)) - (trainNN/2)*log10(2*pi);
    %refine length scale
    lprev = kparams(1);
    kparams(1) = lprev*((rand(1)*0.4)+0.8); %perturb length a little bit
    for oo = 1:100
        %recompute our covariance and likelihood
        K_train_train = GaussianKernel(XTrain,XTrain,kparams);
        Kns = K_train_train + sigman*eye(trainNN);
        plog = -0.5*YTrain'*(Kns\YTrain) - 0.5*log10(det(Kns)) - (trainNN/2)*log10(2*pi);

        %gradient ascent (maximize log likelihood)
        differp = plog - plog_prev;
        differk = kparams(1) - lprev;
        partial = differp / differk; %change in prob / change in params
        lprev = kparams(1);
        kparams(1) = lprev + partial * 0.001;

        %update plog
        plog_prev = plog;
    end
end


%WE STORE THIS FOR ONLINE SHIT
K_inv = inv(Kns);


%store everything in model
Model.K_inv = K_inv;
Model.XTrain = XTrain;
Model.YTrain = YTrain;
Model.kparams = kparams;
Model.sigman = sigman;
Model.dir = [-1,1];
Model.shift = [];
Model.scale = [];
Model.f_classify = @sign; %@round

end