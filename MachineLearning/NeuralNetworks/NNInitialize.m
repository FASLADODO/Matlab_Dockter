function [NN] = NNInitialize(X,Y,layers,nodes)
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
    NN.f_active = @sigmoidFunction;
    NN.f_derivative = @sigmoidDerivative;

    %Initialize NN nerual network weights and bias
    for ll = 1:NN.Layers
       %weight matrix and bias vector for each layer (W_L)
       NN.layer{ll}.weights = rand(NN.Nodes(ll), NN.PrevNodes(ll) );
       NN.layer{ll}.bias = zeros(NN.Nodes(ll),1);
    end

end