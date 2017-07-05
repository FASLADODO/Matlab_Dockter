%Chapter 4 does logistic regression with gradient descent


%Make some data
nn = 1000;

xo = [ones(nn,1),linspace(0,10,nn)'];

params1 = [2;3];
params2 = [2.1;3.5];

y1 = xo*params1 + randn(nn,1)*0.3;
y2 = xo*params2 + randn(nn,1)*0.3;

Data = [xo,y1;xo,y2];
Labels = [ones(nn,1)*0;ones(nn,1)*1];

figure
gscatter(Data(:,2),Data(:,3),Labels)

%% do gradient descent

[NN,SS] = size(Data);

% define logistic
g = @(z) (1 ./ (1 + exp(-z) ) );

theta = randn(SS,1);

cycles = 1000;
thresh = 0.001;
alpha = 0.01; %learning rate

for ii = 1:cycles
    for jj = 1:NN
        X = Data(jj,:);
        P = g(X*theta);
        theta = theta + alpha*(Labels(jj,:) - P)*X';
    end
    E = mean(abs(Labels - g(Data*theta)))
    if(E < thresh)
       break; 
    end
end

theta
finalLabels = g(Data*theta);

figure
plot(1:NN,finalLabels)