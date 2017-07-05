function [outputValue] = NNOnline(NN, inputs)
%online estimate for every row in 'inputs'
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
        layerout{end+1} = NN.f_active( z );
    end
    
    %round output nodes to 0-1
    %layerout{end}
    outputValue(tt,:) = round(layerout{end}); %only the output layers
end

end