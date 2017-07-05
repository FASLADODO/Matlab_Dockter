function [NN,MSE] = NNRegressionTrain(X,Y,NN,alpha,max_epoch,mse_target)
%Based on http://neuralnetworksanddeeplearning.com/chap2.html
% X: Data matrix with rows as samples, columns as dimensions
% Y: Ouput Vector (continuous)
% NN: Neural Network Object (See NNRegressionInitialize.m)
% alpha = 0.1; %learning rate (0 - 1)
% max_epoch = 10000; maximum number of epcohs
% mse_target = 0.01; desired MSE accuracy

%Get sizes
[NNi,SSi] = size(X);
[NNo,SSo] = size(Y);

%default NN structure (1 hidden layers, # nodes = # inputs)
if(isempty(NN))
    [NN] = NNRegressionInitialize(X,Y,1,SSi);
end

%make sure our dimensions are cool
if(SSo ~= NN.Nodes(end))
   error('number of output nodes must equal size of Output data'); 
end
if(NN.Layers ~= length(NN.Nodes))
   error('number of layers must equal number of entries in nodes'); 
end
if(NNo ~= NNi)
    error('#inputs must match #outputs'); 
end

%scale our input data -std<x<std
[Xhat,shift,scale] = MeanVarianceScale(X);
NN.InputScale = scale;
NN.InputShift = shift;

%scale output -std<y<std
[Yhat,shifty,scaley] = NN.f_scale(Y);
NN.OutputShift = shifty; %min(Y); %zeros(1,SSo);
NN.OutputScale = scaley; %ones(1,SSo); %range(Y);

%check if we're changing
RMSE_last = 0;

%iterate backpropogate
STOP_CONDITIONS = false;
epoch_count = 0;
while(~STOP_CONDITIONS)
    
    NNPrev.layer = NN.layer;
    
    %loop through all samples in our training data
    for tt = 1:NNi
        %get current values (we'll use column vectors for math)
        sampleinput = Xhat(tt,:)';
        sampleoutput = Yhat(tt,:)';

        %test the current output with our weights
        layerout = [];
        weightedinput = [];
        %(layerout and weightedinput indices are shifted up one)
        layerout{1} = sampleinput; 
        weightedinput{1} = sampleinput;
        for ll = 1:NN.Layers
            %weighted input at each node in the layer (eq. 25)
            z = NN.layer{ll}.weights * layerout{end} + NN.layer{ll}.bias;
            %activation function is sigmoid for now
            nodeoutput = NN.layer{ll}.f_active(z);

            %store previous layers outputs and sum at each node
            layerout{end+1} = nodeoutput; %output of the current layer (all nodes)
            weightedinput{end+1} = z; % sum value (prior to activation)
        end

        %calculate error at output layer (a_L - y)
        res = layerout{end} - sampleoutput;
        
        %error in the output layer (Eq. BP1)
        delta{NN.Layers} = res .* NN.layer{end}.f_derivative( weightedinput{end} ) ;

        %backpropogate error to previous layers (Eq. BP2)
        %(now we shift our indices back to layer based (1,2,3...)
        for ll = (NN.Layers):-1:2
            %weights in next layer
            w_l1 = NN.layer{ll}.weights;
            %delta in next layer
            delta_l1 = delta{ll};
            %derivative of sigmoid for this layer
            sigz = NN.layer{ll-1}.f_derivative( weightedinput{ll} );
            %updated error delta for this layer
            delta{ll-1} = (w_l1'*delta_l1) .* sigz;
        end
        
        %update bias and weights using error delta terms
        for ll = 1:NN.Layers
            %gradient steps for weights and biases
            %here layerout{ll} is the activations from previous layer
            %because we shifted our indices
            dc_db = delta{ll}; %(Eq. 31)
            dc_dw = delta{ll} * layerout{ll}'; %(Eq. 32)
            
            %actually update them
            NN.layer{ll}.weights = NN.layer{ll}.weights - dc_dw*alpha;
            NN.layer{ll}.bias = NN.layer{ll}.bias - dc_db*alpha;
        end
        
        
    end
    
    
    %test our outputs for this epoch
    Y_Estimate = NNRegressionOnline(NN, X);
    
    %normalized MSE
    ERR = Y - Y_Estimate;
    MSE = mean(NormRowWise(ERR));
    
    %compute change in weights
    for ll = 1:NN.Layers
        diffw = NN.layer{ll}.weights - NNPrev.layer{ll}.weights;
        WeightChange(ll) = mean(mean(abs(diffw)));
    end
    
    if(mean(WeightChange) == 0 && RMSE_last == MSE )
        STOP_CONDITIONS = true;
        fprintf(' Were Not Getting Anywhere\n');
    end
    RMSE_last = MSE;
    
    %increment epochs
    epoch_count = epoch_count + 1;
    
    fprintf('Epoch %d of %d, ', epoch_count,max_epoch);
    fprintf('MSE: %f \n', MSE);
    printf('dW: %f \n', mean(WeightChange));

    if(epoch_count > max_epoch || MSE < mse_target)
        STOP_CONDITIONS = true;
    end
end

fprintf('Total Epochs: %d, ', epoch_count);
fprintf('Final MSE: %f \n', MSE);


end