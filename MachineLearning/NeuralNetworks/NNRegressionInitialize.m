function [NN] = NNRegressionInitialize(X,Y,layers,nodes)
%layers is the number of hidden layers
%nodes is the number of nodes in each hidden layer
%http://neuralnetworksanddeeplearning.com/chap2.html

    
    %get input and output size
    [~,inputnodes] = size(X);
    [~,outputnodes] = size(Y);

    %add in input and output
    layers = layers + 1;
    nodes = [nodes,outputnodes];
    prevnodes = [inputnodes,nodes];

    %stash all our parameters
    NN.Inputs = inputnodes;
    NN.Layers = layers;
    NN.Nodes = nodes;
    NN.PrevNodes = prevnodes;
    NN.InputScale = ones(1,inputnodes);
    NN.InputShift = zeros(1,inputnodes);
    NN.OutputScale = ones(1,outputnodes);
    NN.OutputShift = zeros(1,outputnodes);
    
    %Initialize NN nerual network weights and bias
    for ll = 1:NN.Layers
       %weight matrix and bias vector for each layer (W_L)
       W_scale = 0.1*sqrt(2/NN.PrevNodes(ll)); %to keep weights low when many layers
       NN.layer{ll}.weights = randn(NN.Nodes(ll), NN.PrevNodes(ll) ) * W_scale;
       NN.layer{ll}.bias = zeros(NN.Nodes(ll),1);
       
       %different activations for different layers
       NN.layer{ll}.f_active = @sigmoidFunction;
       NN.layer{ll}.f_derivative = @sigmoidDerivative;
    end
    NN.f_scale = @IdentityScale;
    
    %final layer should be linear
    NN.layer{end}.f_active = @IdentityFunction;
    NN.layer{end}.f_derivative = @IdentityDerivative;
    
    
    %change me to use whatever activation/derivative you want
%     NN.f_active = @SoftPlusFunction;
%     NN.f_derivative = @SoftPlusDerivative;
%     NN.f_scale = @SoftPlusScale;
%     NN.f_active = @sigmoidFunction;
%     NN.f_derivative = @sigmoidDerivative;
%     NN.f_scale = @SigmoidScale;
%     NN.f_active = @BentIdentityFunction;
%     NN.f_derivative = @BentIdentityDerivative;
%     NN.f_scale = @BentIdentityScale;
%     NN.f_active = @TanHFunction;
%     NN.f_derivative = @TanHDerivative;
%     NN.f_scale = @TanHScale;
%     NN.f_active = @ReLUFunction;
%     NN.f_derivative = @ReLUDerivative;
%     NN.f_scale = @ReLUScale;
%     NN.f_active = @LeakyReLUFunction;
%     NN.f_derivative = @LeakyReLUDerivative;
%     NN.f_scale = @LeakyReLUScale;
%     NN.f_active = @IdentityFunction;
%     NN.f_derivative = @IdentityDerivative;
%     NN.f_scale = @IdentityScale;

end