function [network, loss] = NN_Backpropogate(input,output,network)
    %https://sudeepraja.github.io/Neural/
    
    %compute estimates at each layer
    layer_outputs{1} = input;
    layer_inputs{1} = input;
    for i = 1:network.nb_layers
        %raw multiplication
        raw_out = network.layers{i}.weights*layer_outputs{i} + network.layers{i}.biases;
        
        %activation function
        layer_outputs{i+1} = network.layers{i}.activation_fn(raw_out);
        
        %update inputs for next layer
        layer_inputs{i+1} = raw_out;
    end


    %initialize loss at output layer
    delta{network.nb_layers} = network.res_fn(output, layer_outputs{end}) .* network.layers{end}.derivative_fn( layer_inputs{end} ) ;
    
    %now backpropogate the delta terms
    for i = network.nb_layers:-1:2
        %compute delta for previous layer
        delta{i-1} = (network.layers{i}.weights' *delta{i}) .* network.layers{i}.derivative_fn(layer_inputs{i});
    end

    %update weights and biases 
    for i = 1:network.nb_layers
        network.layers{i}.weights = network.layers{i}.weights - network.learning_rate*delta{i}*layer_outputs{i}';
        network.layers{i}.biases = network.layers{i}.biases - network.learning_rate*delta{i};
    end

    %determine new loss
    estimate = NN_Predict(input, network);
    loss = network.loss_fn(output,estimate);
    
end