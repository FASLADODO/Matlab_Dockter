%https://stevenmiller888.github.io/mind-how-to-build-a-neural-network/
%http://neuralnetworksanddeeplearning.com/chap2.html

nn = 100;

d1 = randn(nn,2);
d2 = randn(nn,2) + repmat([2.5,2.5],nn,1);


X = [d1;d2];
Y = [ones(nn,1)*0; ones(nn,1)*1 ];

figure
gscatter(X(:,1),X(:,2),Y)

%% try it with script first

NN = [];

alpha = 0.1; %gradient rate
NN.Inputs = size(X,2);
NN.Layers = 3;
NN.Nodes = [3,5,1]; %3 nodes on hidden layer, 1 output
NN.PrevNodes = [NN.Inputs,NN.Nodes]; %previous layers nodes

%loop through all nodes and layer
for ll = 1:NN.Layers
   %weight matrix and bias vector for each layer (W_L)
   NN.layer{ll}.weights = rand(NN.Nodes(ll), NN.PrevNodes(ll) );
   NN.layer{ll}.bias = zeros(NN.Nodes(ll),1);
end

[NNi,SSi] = size(X);

%we'll mean variance scale the data and sace that scaling for later
[Xhat,shift,scale] = MeanVarianceScale(X);
NN.InputScale = scale;
NN.InputShift = shift;

epochs = 100;
for oo = 1:epochs
    
    %for estimating RMSE and cost
    outputestimate = [];
    costSum = 0;
    
    for tt = 1:NNi
        %get current values (we'll use column vectors for math)
        sampleinput = Xhat(tt,:)';
        sampleoutput = Y(tt,:)';

        %test the current output with our weights
        layerout = [];
        weightedinput = [];
        %(layerout and weightedinput indices are shifted up one)
        layerout{1} = sampleinput; 
        weightedinput{1} = sampleinput;
        for ll = 1:NN.Layers

            %weighted input at each node in the layer (eq. 25)
            z = NN.layer{ll}.weights * layerout{end} + NN.layer{ll}.bias;
            %activation function is sigmoid for now
            nodeoutput = sigmoidFunction(z);

            %store previous layers outputs and sum at each node
            layerout{end+1} = nodeoutput; %output of the current layer (all nodes)
            weightedinput{end+1} = z; % sum value (prior to activation)
        end
        outputestimate(tt,:) = layerout{end}';

        %calculate error at output layer (a_L - y)
        res = layerout{end} - sampleoutput;
        
        %keep track of summation of error cost function (eq. 27)
        costSum = costSum + 0.5*norm(res);
        
        %error in the output layer (Eq. BP1)
        delta{NN.Layers} = res .* sigmoidDerivative( weightedinput{end} ) ;

        %backpropogate error to previous layers (Eq. BP2)
        %(now we shift our indices back to layer based (1,2,3...)
        for ll = (NN.Layers):-1:2
            %weights in next layer
            w_l1 = NN.layer{ll}.weights;
            %delta in next layer
            delta_l1 = delta{ll};
            %derivative of sigmoid for this layer
            sigz = sigmoidDerivative( weightedinput{ll} );
            %updated error delta for this layer
            delta{ll-1} = (w_l1'*delta_l1) .* sigz;
        end
        
        %update bias and weights using error delta terms
        for ll = 1:NN.Layers
            %gradient steps for weights and biases
            %here layerout{ll} is the activations from previous layer
            %because we shifted our indices
            dc_db = delta{ll}; %(Eq. 31)
            dc_dw = delta{ll} * layerout{ll}'; %(Eq. 32)
            
            %actually update them
            NN.layer{ll}.weights = NN.layer{ll}.weights - dc_dw*alpha;
            NN.layer{ll}.bias = NN.layer{ll}.bias - dc_db*alpha;
        end
    end

    %compute errors / cost
    C = (1/NNi) * costSum;
    rmse = sqrt(mean( (outputestimate - Y).^2 ));

end

Y_est = round(outputestimate);

figure
gscatter(X(:,1),X(:,2),Y_est)


%% Now we weaponized it

layers = 1;
nodes = 3;
alpha = 0.5;

clear NN
[NN] = NNInitialize(X,Y,layers,nodes);

[NN,MSE] = NNTrain(X,Y,NN,alpha,100,0.01);

Yest = NNOnline(NN, X);
    
figure
gscatter(X(:,1),X(:,2),Yest)

