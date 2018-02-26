function [layer] = NN_Layer(inputs, outputs, activation)
    %define sizes
    layer.inputs = inputs;
    layer.outputs = outputs;
    layer.activation_str = activation;
    layer.scale = 0.001;

    %initialize weights
    layer.weights = (rand(outputs,inputs)-0.5)*2.0*layer.scale;
    
    %initialize biases (rand, zeros)
    layer.biases = (rand(outputs,1)-0.5)*2.0*layer.scale;
    
    %activation functions
    if( strcmp(activation,'relu') == 1)
       %relu activation
       layer.activation_fn =  @(x) max(x,zeros(size(x)) );
       layer.derivative_fn = @(x) double(x > 0);
    elseif( strcmp(activation,'sigmoid') == 1)
       %sigmoid activation
       layer.activation_fn =  @(x) 1.0 ./ (1.0 + exp(-x));
       layer.derivative_fn = @(x) (1.0 ./ (1.0 + exp(-x))) .* (1.0 - (1.0 ./ (1.0 + exp(-x))) );
    else
       %default is linear
       layer.activation_str = 'linear';
       layer.activation_fn =  @(x) x;
       layer.derivative_fn = @(x) 1.0;
    end
end
