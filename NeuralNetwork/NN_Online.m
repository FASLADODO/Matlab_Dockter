function outputNode = NN_Online(NN, inputs)
%Rods attempt at a neural network

layers = length(NN.layer);

nodes.layer{1} = inputs;
for ll = 1:layers
    for jj = 1:length(NN.layer{ll}.node)
        %t = bsxfun(@times, params.weights{ll}, node_in)
        nodes.layer{ll+1}(jj) = sigmoidal( sum(NN.layer{ll}.node{jj}.weights .* nodes.layer{ll}) - NN.layer{ll}.node{jj}.threshold );
    end
end

outputNode = nodes.layer;

end