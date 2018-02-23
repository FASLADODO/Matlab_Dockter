function [output] = NN_Predict(inputs, network)
    %forward prediction through network
    for i = 1:network.nb_layers
        %raw multiplication
        raw_out = network.layers{i}.weights*inputs + network.layers{i}.biases;
        
        %activation function
        output = network.layers{i}.activation_fn(raw_out);
        
        %update inputs for next layer
        inputs = output;
    end
end