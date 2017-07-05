function [outputValue] = NNRegressionOnline(NN, inputs)
%online regression output for every row in 'inputs'
%uses model trained 'NN'
%http://neuralnetworksanddeeplearning.com/chap2.html

%scale and shift inputs
[NNi,~] = size(inputs);
inputs = (inputs - repmat(NN.InputShift,NNi,1) ) ./  repmat(NN.InputScale,NNi,1);

outputValue = [];
for tt = 1:NNi
    layerout = [];
    layerout{1} = inputs(tt,:)';
    for ll = 1:NN.Layers
        %weighted input at each node in the layer (eq. 23)
        z = NN.layer{ll}.weights * layerout{end} + NN.layer{ll}.bias;
        %output activation function (eq. 25)
        layerout{end+1} = NN.layer{ll}.f_active( z );
    end
    
    %scale output back to regression values
    %Order of operations matters!
    outputValue(tt,:) = (layerout{end} .* NN.OutputScale ) + NN.OutputShift; %only the output layers
end

end