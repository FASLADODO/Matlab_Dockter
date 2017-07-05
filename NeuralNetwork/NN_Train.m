function NN = NN_Train(Training) 
    %Taken from here
    %https://takinginitiative.wordpress.com/2008/04/03/basic-neural-network-tutorial-theory/
    
    %VARIABLES (todo)
    alpha = 0.3; %learning rate (0 - 1)
    max_epoch = 10000;
    accuracy_target = 0.005;

    %initialize neural net
    training_set = length(Training);
    layers = 2; %input, hidden, output
    hidden_nodes = 8; %can be changed
    input_nodes = length(Training{1}.inputs);
    output_nodes = length(Training{1}.outputs);
    nodes = [ hidden_nodes, output_nodes];
    weight_size = [ input_nodes, hidden_nodes];
    scale  = 2.4 / input_nodes;
    
    
    %create initial struct
    NN = NN_Create(layers,nodes);
    

    %Initialize random numbers for weights
    for ll = 1:layers
        for jj = 1:nodes(ll)
            NN.layer{ll}.node{jj}.weights = 2 *scale * rand(1,weight_size(ll)) - scale;
            NN.layer{ll}.node{jj}.threshold = 1; %scale * rand(1);
        end
    end
    
    %iterative approach to get proper weights (gradient descent)
    STOP_CONDITIONS = true;
    epoch_count = 0;
    while(STOP_CONDITIONS)
        %For MSE
        MSE_numerator = 0;
        
        %output error gradients
        for tt = 1:training_set
            
            %Feed forward
            outputNode = NN_Online(NN, Training{tt}.inputs);
            
            %Zero out
            deltaHiddenOutput = zeros(hidden_nodes,output_nodes);
            deltaInputHidden = zeros(input_nodes,hidden_nodes);
            
            %get output layer gradients
            for k = 1:output_nodes
                outputErrorGradients(k) = outputNode{end}(k)*(1 - outputNode{end}(k)).*(Training{tt}.outputs(k) - outputNode{end}(k));

                for j = 1:hidden_nodes
                    deltaHiddenOutput(j,k) = alpha * outputNode{end-1}(j) * outputErrorGradients(k);
                end
                
            end
            
            %hidden error gradients
            for j = 1:hidden_nodes
                weightedSum = 0;
                for k = 1:output_nodes
                    % sum up: w_{jk}*delta_{k}
                    weightedSum =  weightedSum + NN.layer{end}.node{k}.weights(j) * outputErrorGradients(k);
                end
                %delta_j = yj(1-yj)sum(wjk*deltak)
                hiddenErrorGradients(j) = outputNode{end-1}(j) * ( 1 - outputNode{end-1}(j) ) * weightedSum;
                
                for i = 1:input_nodes
                    %calculate change in weight 
                    deltaInputHidden(i,j) = deltaInputHidden(i,j) + alpha * outputNode{1}(i) * hiddenErrorGradients(j); 
                end
            end
            
            %Update the weights 
            for i = 1:input_nodes
                for j = 1:hidden_nodes
                    NN.layer{1}.node{j}.weights(i) = NN.layer{1}.node{j}.weights(i) + deltaInputHidden(i,j);
                end
            end
            for j = 1:hidden_nodes
                for k = 1:output_nodes
                    NN.layer{2}.node{k}.weights(j) = NN.layer{2}.node{k}.weights(j) + deltaHiddenOutput(j,k);
                end
            end
            
            %test accuracy
            outputtest = NN_Online(NN, Training{tt}.inputs);
            errortest = norm(outputtest{end} - Training{tt}.outputs);
            MSE_numerator = MSE_numerator + errortest;
        end
        %normalize MSE
        MSE = MSE_numerator / training_set;
        %increment epochs
        epoch_count = epoch_count + 1;
        
        if(epoch_count > max_epoch || MSE < accuracy_target)
            STOP_CONDITIONS = false;
        end
    end
    epoch_count
end