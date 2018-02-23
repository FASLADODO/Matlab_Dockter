function [network] = NN_Compile(loss, learning_rate, varargin)
    %network parameters
    network.nb_layers = nargin - 2;
    network.learning_rate = learning_rate;
    network.loss_str = loss;
    
     if( strcmp(loss,'cross_entropy') == 1)
       %cross entropy loss (target, estimate)
       network.res_fn =  @(t,y) -t .* log(y);
       network.loss_fn =  @(t,y) sum( -t .* log(y)) / numel(y);
     elseif( strcmp(loss,'mse') == 1)
       %mean sqaured error  (target, estimate)
       network.res_fn =  @(t,y) (t-y);
       network.loss_fn =  @(t,y) sum(t-y) / numel(y);
     end
    
    %network layers
    for i = 1:network.nb_layers
        network.layers{i} = varargin{i};
    end
end