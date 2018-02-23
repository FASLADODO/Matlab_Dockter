%example neural network

inputs = 2;
hiddensize1 = 6;
outputs = 2;
learning_rate = 0.01;
nb_epochs = 10000;

%define layers
layer1 = NN_Layer(inputs, hiddensize1, 'relu');
layer2 = NN_Layer(hiddensize1, outputs, 'sigmoid');


%define outputs
network = NN_Compile('cross_entropy', learning_rate, layer1, layer2);

%sample inputs xor problem
nsamples = 4;
x_data = [0, 0, 1, 1;
          0, 1, 0, 1];
y_data = [0, 1, 1, 0];

for i = 1:nb_epochs
    
    loss = 0;
    for j = 1:nsamples
        sample = x_data(:,j);
        label = y_data(j);
    
        [network, temploss] = NN_Backpropogate(sample,label,network);
        loss = loss + temploss;
    end
    fprintf('Epoch %d of %d, ', i, nb_epochs);
    fprintf('Current Loss: %f \n', loss / nb_epochs);
end


%do forward prediction
estimate = NN_Predict([1;0], network)
